#include "parsebundler.h"

#include "opencv2/opencv.hpp"

#include <string.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <utility>
#include <numeric>
#include <algorithm>

//using namespace std;
using std::string;
using std::ifstream;
using std::ofstream;
using std::istringstream;
//using std::cout;
using std::endl;


PARSE_BUNDLER::PARSE_BUNDLER()
{
	mNumbPoints = 0;
	mNumCameras = 0;
}

PARSE_BUNDLER::~PARSE_BUNDLER()
{
	mAll_pic_cameras.ClearPicsCameras();
	mFeature_infos.clear();
}


size_t PARSE_BUNDLER::GetNumPoints() const{
	return mNumbPoints;
}

size_t PARSE_BUNDLER::GetNumCameras() const{
	return mNumCameras;
}


//clear all the data
void  PARSE_BUNDLER::ClearData(){
	mAll_pic_cameras.ClearPicsCameras();
	mFeature_infos.clear();
	mNumbPoints = 0;
	mNumCameras = 0;
}


//parse the bundler file
/*Each camera entry contains the estimated camera intrinsics and extrinsics,
  and has the form:
<f> <k1> <k2>   [the focal length, followed by two radial distortion coeffs]
<R>             [a 3x3 matrix representing the camera rotation]
<t>             [a 3-vector describing the camera translation]
Each point entry has the form
<position>      [a 3-vector describing the 3D position of the point]
<color>         [a 3-vector describing the RGB color of the point]
<view list>     [a list of views the point is visible in]
The view list begins with the length of the list 
(i.e., the number of cameras the point is visible in).
The list is then given as a list of quadruplets <camera> <key> <x> <y>
The pixel positions are floating point numbers in a coordinate system 
where the origin is the center of the image, the x-axis increases 
to the right, and the y-axis increases towards the top of the image. 
Thus, (-w/2, -h/2) is the lower-left corner of the image, 
and (w/2, h/2) is the top-right corner (where w and h are the width 
and height of the image).
*/

bool PARSE_BUNDLER::ParseBundlerFile( const std::string& bundle_file )
{
  std::ifstream instream( bundle_file, std::ios::in );
	if (!instream.is_open()){
    std::cout << "open bundler fail: " << bundle_file << std::endl;
		return 0;
	}

	string line_in_file;
	getline(instream, line_in_file);//header

	instream >> mNumCameras >> mNumbPoints;

	auto & mCameras = mAll_pic_cameras.GetAllCameras();
	//mCameras.clear();
	mCameras.resize(mNumCameras);
	//mFeature_infos.clear();
	mFeature_infos.resize(mNumbPoints);

	//load the camera parameters 
	for (int i = 0; i < mNumCameras; i++)
	{
		instream >> mCameras[i].focal_length 
			>> mCameras[i].k1 >> mCameras[i].k2;
		
		instream >> mCameras[i].rotation(0, 0)
			>> mCameras[i].rotation(0, 1)
			>> mCameras[i].rotation(0, 2)
			>> mCameras[i].rotation(1, 0)
			>> mCameras[i].rotation(1, 1)
			>> mCameras[i].rotation(1, 2)
			>> mCameras[i].rotation(2, 0)
			>> mCameras[i].rotation(2, 1)
			>> mCameras[i].rotation(2, 2);
		
		instream >> mCameras[i].translation(0)
			>> mCameras[i].translation(1)
			>> mCameras[i].translation(2);
		//mCameras[i].id = i;
	}

	//load the points
	int r, g, b;
	for (int i = 0; i < mNumbPoints; i++)
	{
		instream >> mFeature_infos[i].mPoint.x
			>> mFeature_infos[i].mPoint.y
			>> mFeature_infos[i].mPoint.z
			>> r >> g >> b;
		mFeature_infos[i].mPoint.r = (unsigned char)r;
    mFeature_infos[i].mPoint.g = (unsigned char)g;
    mFeature_infos[i].mPoint.b = (unsigned char)b;
		
		int view_lenth = 0;
		instream >> view_lenth;
		mFeature_infos[i].mView_list.resize(view_lenth);

		for (int j = 0; j < view_lenth; j++)
		{
			instream >> mFeature_infos[i].mView_list[j].camera
				>> mFeature_infos[i].mView_list[j].key
				>> mFeature_infos[i].mView_list[j].x
				>> mFeature_infos[i].mView_list[j].y;
		}
	}

	instream.close();
	return 1;
}


//convert the bundle into ply file
static char ply_header[] =
"ply\n"
"format ascii 1.0\n"
"element vertex %d\n"
"property float x\n"
"property float y\n"
"property float z\n"
"property uchar diffuse_red\n"
"property uchar diffuse_green\n"
"property uchar diffuse_blue\n"
"end_header\n";
bool PARSE_BUNDLER::ConvertPly( const std::string& ply_file ) const
{
  FILE *f = fopen( ply_file.c_str(), "wt" );
  if ( f == NULL ) {
    printf( "Error opening file %s for writing\n", ply_file.c_str() );
    return 0;
  }
  /* save only 3d points with more than 2 feature points*/
  size_t cnt = 0;
  for ( const auto& point_info : mFeature_infos ){
    if ( point_info.mView_list.size() > 2 ){
      cnt += point_info.mView_list.size();
    }
  }

  /* Print the ply header */
  fprintf( f, ply_header, cnt );
  for ( const auto& point_info : mFeature_infos ){
    if ( point_info.mView_list.size( ) > 2 )
    fprintf( f, "%0.6e %0.6e %0.6e %d %d %d\n",
      point_info.mPoint.x, point_info.mPoint.y, point_info.mPoint.z,
      point_info.mPoint.r, point_info.mPoint.g, point_info.mPoint.b );
  }

  fclose( f );
  return 1;
}
/*
//convert the bundle into ply file
bool PARSE_BUNDLER::ConvertPly( const std::string& ply_file ) const{
  std::ofstream outfile( ply_file, std::ios::out  | std::ios::trunc );
  if ( !outfile.is_open( ) ){
    std::cout << "open bundler fail: " << ply_file << std::endl;
    return 0;
  }
  outfile << "ply"                            << std::endl;
  outfile << "format ascii 1.0"               << std::endl;
  outfile << "element vertex " << mNumbPoints << std::endl;
  outfile << "property float x"               << std::endl;
  outfile << "property float y"               << std::endl;
  outfile << "property float z"               << std::endl;
  outfile << "property uchar diffuse_red"     << std::endl;
  outfile << "property uchar diffuse_green"   << std::endl;
  outfile << "property uchar diffuse_bule"    << std::endl;
  outfile << "end_header"                     << std::endl;
  for (const auto& point_info : mFeature_infos ){
    outfile << point_info.mPoint.x << " "
      << point_info.mPoint.y << " "
      << point_info.mPoint.z << " "
      << (int)point_info.mPoint.r << " "
      << (int)point_info.mPoint.g << " "
      << (int)point_info.mPoint.b << " "
      << std::endl;
  }
  return 1;
}*/


// do some statistics of 3d points
// the length of view list and point color
// format: mNumCameras mNumbPoints
// format: length rgb_bundler  mean_rgb var_rgb
void PARSE_BUNDLER::DoStatPoints( const std::string& stat_res_file) const
{
  std::ofstream outfile( stat_res_file, std::ios::trunc );
  if ( !outfile.is_open() ){
    std::cout << "Open file fail: " << stat_res_file << endl;
    return;
  }

  outfile << mNumCameras << " " << mNumbPoints << std::endl;

  size_t total_view_list_len = 0;
  std::vector<unsigned int> r,g,b;

  for ( const auto& point_info : mFeature_infos )
  {
    r.clear(); g.clear(); b.clear();
    total_view_list_len += point_info.mView_list.size();
    
    for ( const auto& view : point_info.mView_list )
    {
      r.push_back( view.rgb[0] );
      g.push_back( view.rgb[1] );
      b.push_back( view.rgb[2] );
    }
    double rgb_mean[3];
    double rgb_variance[3];
    rgb_mean[0] = std::accumulate( std::begin( r ), std::end( r ), 0.0 );
    rgb_mean[1] = std::accumulate( std::begin( g ), std::end( g ), 0.0 );
    rgb_mean[2] = std::accumulate( std::begin( b ), std::end( b ), 0.0 );
    rgb_mean[0] = rgb_mean[0] / r.size();
    rgb_mean[1] = rgb_mean[1] / g.size( );
    rgb_mean[2] = rgb_mean[2] / b.size( );

    rgb_variance[0] = rgb_variance[1] = rgb_variance[2] = 0;
    std::for_each( std::begin( r ), std::end( r ), [&]( const float d ){
      rgb_variance[0] += ( d - rgb_mean[0] )*( d - rgb_mean[0] );
    } );
    std::for_each( std::begin( g ), std::end( g ), [&]( const float d ){
      rgb_variance[1] += ( d - rgb_mean[1] )*( d - rgb_mean[1] );
    } );
    std::for_each( std::begin( b ), std::end( b ), [&]( const float d ){
      rgb_variance[2] += ( d - rgb_mean[2] )*( d - rgb_mean[2] );
    } );
    rgb_variance[0] = std::sqrt( rgb_variance[0] / r.size() );
    rgb_variance[1] = std::sqrt( rgb_variance[1] / g.size() );
    rgb_variance[2] = std::sqrt( rgb_variance[2] / b.size() );

    outfile << point_info.mView_list.size( ) << " "
      << point_info.mPoint.x << " "
      << point_info.mPoint.y << " "
      << point_info.mPoint.z << " "
      << (int)point_info.mPoint.r << " "
      << (int)point_info.mPoint.g << " "
      << (int)point_info.mPoint.b << " "
      << (int)rgb_mean[0] << " "
      << (int)rgb_mean[1] << " "
      << (int)rgb_mean[2] << " "
      << (int)rgb_variance[0] << " "
      << (int)rgb_variance[1] << " "
      << (int)rgb_variance[2] << " "
      << std::endl;
  }
  std::cout << "total 3d points: "  << mNumbPoints << std::endl;
  std::cout << "total cameras: "    << mNumCameras << std::endl;
  std::cout << "avg cameras each 3d points in: " << total_view_list_len / mNumbPoints << std::endl;

}

// load feature point color from pictures
bool PARSE_BUNDLER::LoadFeatPointRGB(){
  auto& picture = mAll_pic_cameras.GetAllPictures();
#pragma omp parallel for
  for ( int i = 0; i < mNumbPoints; ++i )
  {
    for ( int j = 0; j < mFeature_infos[i].mView_list.size(); ++j )
    {
      auto & view = mFeature_infos[i].mView_list[j];
      auto & keypoint_vec = picture[view.camera].GetFeaturePoint( );

      //use the SIFT x, y, of which the origin is the left-up corner
      view.rgb[0] = keypoint_vec[view.key].rgb[0];
      view.rgb[1] = keypoint_vec[view.key].rgb[1];
      view.rgb[2] = keypoint_vec[view.key].rgb[2];
    }
  }
  return 1;
}

//load the .key info 
bool PARSE_BUNDLER::LoadCameraInfo( bool loadDesc )
{
	auto& picture = mAll_pic_cameras.GetAllPictures();
	//make sure that database #pictures equal to #cameras
	assert(picture.size() == mAll_pic_cameras.GetAllCameras().size());

#pragma omp parallel for
	for (int i = 0; i < mNumbPoints; ++i)
	{
		mFeature_infos[i].mDescriptor.resize(mFeature_infos[i].mView_list.size());
		for (int j = 0; j < mFeature_infos[i].mView_list.size(); ++j)
		{
			auto & view = mFeature_infos[i].mView_list[j];
			auto & keypoint_vec = picture[mFeature_infos[i].mView_list[j].camera].GetFeaturePoint();

			//use SIFT x, y coordinates, of which the origin is the left-up corner
      view.x            = keypoint_vec[view.key].x;
      view.y            = keypoint_vec[view.key].y;
			view.scale			  = keypoint_vec[view.key].scale;
			view.orientation	= keypoint_vec[view.key].orientation;
      
      if ( !loadDesc ) continue;
			//mDescriptor: the jth descriptor
			auto & descriptor = picture[mFeature_infos[i].mView_list[j].camera].GetDescriptor();
			mFeature_infos[i].mDescriptor[j] = std::move(descriptor[view.key]);
		}
	}

	return 1;
}


//after load original bundler file, mask the query image
void PARSE_BUNDLER::FindQueryPicture(const std::string& s)
{
	ifstream is(s, std::ios::in);
	if (is.is_open() == 0){
		std::cout << "query_picture_list.txt open fail: "<< s << endl;
		return;
	}

	mPic_query_mask.clear();
  mPic_query_mask.reserve( 1000 );
	string line;
	while (getline(is, line))
	{
		//stringstream istrstream;
		mPic_query_mask.push_back( line[0] == 'q' );
	}

}

// write query image into query bundler
// write point mode: 
//                  0, without point; 
//                  1, point referred in db index; 
//                  2, point referred in query index;
//                  
// mode =2 , means it is independent bundle.out file with small partition of original images
void PARSE_BUNDLER::WriteQueryBundler( const std::string& sBundleQuery, int iWritePointMode, bool bWithPointIndex ) const
{
  std::map<size_t, size_t> mapImgOrginToQuery;
  //std::map<size_t, size_t> mapImgOrginToDB;

  // build the index from the original img index to the query img index
  size_t cntImgQuery = 0;
  //size_t cntImgDB = 0;
  for ( size_t cntImage = 0; cntImage < mNumCameras; cntImage++ )
  {
    if ( true == mPic_query_mask[cntImage] ){
      mapImgOrginToQuery[cntImage] = cntImgQuery++;
    }
    // else{ mapImgOrginToDB[cntImage] = cntImgDB++; }
  }
  std::cout << "num of query image: " << cntImgQuery  << std::endl;
  //std::cout << "num of DB image: "    << cntImgDB     << std::endl;

  // find 3d points remain in bundle.db.out and those only exist in query bundler
  size_t cntPointsInDB = 0;
  size_t cntPointsInQuery = 0;
  std::vector<bool> bPointInDB( mNumbPoints, false );
  std::vector<bool> bPointInQuery( mNumbPoints, false );
  std::map<size_t, size_t> map3DPointsOrigToDB;
  std::map<size_t, size_t> map3DPointsOrigToQuery;

  for ( size_t i = 0; i < mNumbPoints; i++ )
  {
    size_t cntPointHaveDBCams = 0;
    size_t cntPointHaveQueryCames = 0;
    // find if this point should be remained in db
    for ( const auto& view : mFeature_infos[i].mView_list )
    {
      if ( false == mPic_query_mask[view.camera] )
      {
        ++cntPointHaveDBCams;
      }
      else
      {
        ++cntPointHaveQueryCames;
      }

      if ( cntPointHaveDBCams >= 2 && cntPointHaveQueryCames >= 2) break;
    }

    // build point index from original to db
    if ( cntPointHaveDBCams >= 2 )
    {
      bPointInDB[i] = true; // this point will remain in db
      map3DPointsOrigToDB[i] = cntPointsInDB++;
    }
    // build point index from original to query
    if( cntPointHaveQueryCames >=2 )
    { 
      bPointInQuery[i] = true; // this point will remain in query
      map3DPointsOrigToQuery[i] = cntPointsInQuery++;
    }
  }
  std::cout << "num of db points: " << cntPointsInDB << std::endl;
  std::cout << "num of query points: " << cntPointsInQuery << std::endl;

  std::ofstream os( sBundleQuery, std::ios::out | std::ios::trunc );
  if ( false == os.is_open() ){
    std::cout << " open query bundle file fail: " << sBundleQuery << endl;
    return;
  }

  const auto & mCameras = mAll_pic_cameras.GetAllCameras();

  os << "# Bundle file v0.3" << endl;
  size_t numPoint[3] = { 0, cntPointsInDB, cntPointsInQuery };
  size_t cntPoint = numPoint[iWritePointMode];
  os << cntImgQuery << " " << cntPoint << endl;
  os << std::resetiosflags( std::ios::fixed )
    << std::setiosflags( std::ios::scientific )
    << std::setprecision( 9 );
  for ( size_t i = 0; i < mNumCameras; i++ )
  {
    if ( true == mPic_query_mask[i] )
    {
      os << mCameras[i].focal_length << " "
        << mCameras[i].k1 << " " << mCameras[i].k2 << endl;
      os << mCameras[i].rotation( 0, 0 ) << " "
        << mCameras[i].rotation( 0, 1 ) << " "
        << mCameras[i].rotation( 0, 2 ) << endl;
      os << mCameras[i].rotation( 1, 0 ) << " "
        << mCameras[i].rotation( 1, 1 ) << " "
        << mCameras[i].rotation( 1, 2 ) << endl;
      os << mCameras[i].rotation( 2, 0 ) << " "
        << mCameras[i].rotation( 2, 1 ) << " "
        << mCameras[i].rotation( 2, 2 ) << endl;
      os << mCameras[i].translation( 0 ) << " "
        << mCameras[i].translation( 1 ) << " "
        << mCameras[i].translation( 2 ) << endl;
    }
  }

  // without saving points
  if ( 0 == iWritePointMode )
  { 
    os.close();
    return; 
  }

  // save the point and view list
  for ( size_t i = 0; i < mNumbPoints; i++ )
  {
    if ( (1 == iWritePointMode && true == bPointInDB[i])
      || (2 == iWritePointMode && true == bPointInQuery[i] ) )
    {
      size_t view_lenth = mFeature_infos[i].mView_list.size();
      size_t view_lenth_true = 0;

      //calculate true view_length
      for ( size_t j = 0; j < view_lenth; j++ )
      {
        if ( true == mPic_query_mask[mFeature_infos[i].mView_list[j].camera] )
        {
          view_lenth_true++;
        }
      }
      // there is no query cameras then do not save this point
      if ( 0 == view_lenth_true ) continue;

      os << std::resetiosflags( std::ios::fixed )
        << std::setiosflags( std::ios::scientific )
        << std::setprecision( 9 );
      os << mFeature_infos[i].mPoint.x << " "
        << mFeature_infos[i].mPoint.y << " "
        << mFeature_infos[i].mPoint.z << std::endl
        << (int)mFeature_infos[i].mPoint.r << " "
        << (int)mFeature_infos[i].mPoint.g << " "
        << (int)mFeature_infos[i].mPoint.b;

      //if write the index of 3d point
      if ( bWithPointIndex && ( 1 == iWritePointMode ) )
      {
        os << " " << map3DPointsOrigToDB[i];
      }
      else if ( bWithPointIndex && ( 2 == iWritePointMode ) ){ 
        os << " " << map3DPointsOrigToQuery[i];
      }
      os << std::endl << view_lenth_true << " ";

      for ( size_t j = 0; j < view_lenth; j++ )
      {
        os << std::resetiosflags( std::ios::scientific )
          << std::setiosflags( std::ios::fixed )
          << std::setprecision( 4 );
        size_t cntCamera = mFeature_infos[i].mView_list[j].camera;
        if ( true == mPic_query_mask[cntCamera] )
        {
          const auto mapQuery_iter = mapImgOrginToQuery.find( cntCamera );
          if ( mapQuery_iter != mapImgOrginToQuery.cend() )
          { 
            os << mapQuery_iter->second << " "
              << mFeature_infos[i].mView_list[j].key << " "
              << mFeature_infos[i].mView_list[j].x << " "
              << mFeature_infos[i].mView_list[j].y << " ";
          }
          else{
            std::cout << " save view list error. parse_bundler.cpp line 528" << std::endl;
            return;
          }
        }
      }
      os << std::endl;
    }
  }
  os.close();

  /* save the 3d point index map */
  // then save the match information
  std::ofstream ofs( sBundleQuery + ".points.map", std::ios::trunc );
  if ( !ofs.is_open() ){
    std::cout << "Open points.map fail: " << sBundleQuery + ".points.map" << std::endl;
    return;
  }
  for ( auto iter = map3DPointsOrigToDB.cbegin(); iter != map3DPointsOrigToDB.cend(); ++iter )
  {
    ofs << iter->first << " " << iter->second << std::endl;
  }
  ofs.close();

}

// write db image into db bundler
void PARSE_BUNDLER::WriteDBBundler( const std::string& sBundleDB ) const
{
  //std::map<size_t, size_t> mapImgOrginToQuery;
  std::map<size_t, size_t> mapImgOrginToDB;

  // build the index from the original img index to the query img index
  size_t cntImgDB = 0;
  for ( size_t cntImage = 0; cntImage < mNumCameras; cntImage++ )
  {
    if ( false == mPic_query_mask[cntImage] ){
      mapImgOrginToDB[cntImage] = cntImgDB++;
    }
  }
  std::cout << "num of db image: " << cntImgDB << std::endl;

  // find 3d points remain in bundle.db.out
  size_t cntPointsInDB = 0;
  std::vector<bool> bPointInDB( mNumbPoints, false );
  std::map<size_t, size_t> map3DPointsOrigToDB;

  for ( size_t i = 0; i < mNumbPoints; i++ )
  {
    int cntPointHaveDBCams = 0;
    // find if this point should be remained in db
    for ( const auto& view : mFeature_infos[i].mView_list )
    {
      if ( false == mPic_query_mask[view.camera] ){
        ++cntPointHaveDBCams;
        if ( cntPointHaveDBCams >= 2 ) break;
      }
    }

    // build point index from original to db
    if ( cntPointHaveDBCams >= 2 )
    {
      bPointInDB[i] = true; // this point will remain in db
      map3DPointsOrigToDB[i] = cntPointsInDB++;
    }
  }
  std::cout << "num of db points: " << cntPointsInDB << std::endl;

  // then save the bundle.db file
  ofstream os( sBundleDB, std::ios::out | std::ios::trunc );
  if ( false == os.is_open() ){
    std::cout << " open query bundle file fail: " << sBundleDB << endl;
    return;
  }

  const auto & mCameras = mAll_pic_cameras.GetAllCameras();

  os << "# Bundle file v0.3" << endl;

  os << cntImgDB << " " << cntPointsInDB << endl;
  os << std::resetiosflags( std::ios::fixed ) 
    << std::setiosflags( std::ios::scientific ) 
    << std::setprecision( 9 );
  //save the camera information
  for ( size_t i = 0; i < mPic_query_mask.size(); i++ ){
    if ( false == mPic_query_mask[i] ){
      os << mCameras[i].focal_length << " "
        << mCameras[i].k1 << " " << mCameras[i].k2 << endl;
      os << mCameras[i].rotation( 0, 0 ) << " "
        << mCameras[i].rotation( 0, 1 ) << " "
        << mCameras[i].rotation( 0, 2 ) << endl;
      os << mCameras[i].rotation( 1, 0 ) << " "
        << mCameras[i].rotation( 1, 1 ) << " "
        << mCameras[i].rotation( 1, 2 ) << endl;
      os << mCameras[i].rotation( 2, 0 ) << " "
        << mCameras[i].rotation( 2, 1 ) << " "
        << mCameras[i].rotation( 2, 2 ) << endl;
      os << mCameras[i].translation( 0 ) << " "
        << mCameras[i].translation( 1 ) << " "
        << mCameras[i].translation( 2 ) << endl;
    }
  }

  // save the point and view list
  for ( size_t i = 0;  i < mNumbPoints; i++ )
  {
    if ( true == bPointInDB[i] )
    {
      size_t view_lenth = mFeature_infos[i].mView_list.size();
      size_t view_lenth_true = 0;

      //calculate true view_length include only db points in cameras
      for ( size_t j = 0; j < view_lenth; j++ )
      {
        if ( false == mPic_query_mask[mFeature_infos[i].mView_list[j].camera] )
        {
          view_lenth_true++;
        }
      }

      os << std::resetiosflags( std::ios::fixed ) 
        << std::setiosflags( std::ios::scientific ) 
        << std::setprecision( 9 );
      os << mFeature_infos[i].mPoint.x << " "
        << mFeature_infos[i].mPoint.y << " "
        << mFeature_infos[i].mPoint.z << std::endl
        << (int)mFeature_infos[i].mPoint.r << " "
        << (int)mFeature_infos[i].mPoint.g << " "
        << (int)mFeature_infos[i].mPoint.b << std::endl
        << view_lenth_true << " ";

      for ( size_t j = 0; j < view_lenth; j++ )
      {
        os << std::resetiosflags( std::ios::scientific ) 
          << std::setiosflags( std::ios::fixed ) 
          << std::setprecision( 4 );
        if ( false == mPic_query_mask[mFeature_infos[i].mView_list[j].camera] )
        {
          const auto map_iter = mapImgOrginToDB.find( mFeature_infos[i].mView_list[j].camera );
          if ( map_iter != mapImgOrginToDB.cend() )
          { 
            os <<  map_iter->second  << " "
              << mFeature_infos[i].mView_list[j].key << " "
              << mFeature_infos[i].mView_list[j].x << " "
              << mFeature_infos[i].mView_list[j].y << " ";
          }
          else{ 
            std::cout << " save view list err. parse_bundler.cpp line 672" << std::endl;
            return;
          }
        }
      }
      os << std::endl;
    }
  }

  os.close();

  /* save the 3d point index map */
  // then save the match information
  std::ofstream ofs( sBundleDB + ".points.map", std::ios::trunc );
  if ( !ofs.is_open() ){
    std::cout << "Open points.map fail: " << sBundleDB + ".points.map" << std::endl;
    return;
  }
  for (auto iter = map3DPointsOrigToDB.cbegin(); iter != map3DPointsOrigToDB.cend(); ++iter )
  {
    ofs << iter->first << " " << iter->second << std::endl;
  }
  ofs.close();
}


// find the query image true match, work for original bundle.out
void PARSE_BUNDLER::FindQueryFeatureTrueMatch( const std::string& sTrueMatchFile ) const
{
  std::map<size_t, size_t> mapImgOrginToQuery;
  // std::map<size_t, size_t> mapOrginToDB;
  
  // build the index from the original img index to the query img index
  size_t cntImgQuery = 0;
  for ( size_t cntImage = 0; cntImage < mNumCameras; cntImage++ )
  {
    if ( true == mPic_query_mask[cntImage] ){
      mapImgOrginToQuery[cntImage] = cntImgQuery++;
    }
  }
  std::cout << "num of query image: " << cntImgQuery << std::endl;

  // find 3d points only see in query pictures
  size_t cntPointsInDB = 0;
  std::map<size_t, size_t> map3DPointsOrigToDB;
  std::vector<std::vector<std::pair<size_t, size_t>>> vPairImgFeatTo3DPointMatch( cntImgQuery );
  for ( size_t i = 0; i < mNumbPoints; i++ )
  {
    int cntPointHaveDBCams = 0;
    // find if this point should be remained in db
    for ( const auto& view : mFeature_infos[i].mView_list )
    {
      if ( false == mPic_query_mask[view.camera] ){
        ++cntPointHaveDBCams;
        if ( cntPointHaveDBCams >= 2 ) break;
      }
    }
    
    // build point index from original to db
    if ( cntPointHaveDBCams >= 2 )
    {
      map3DPointsOrigToDB[i] = cntPointsInDB;
      size_t kCamQuery = 0;
      // build the query match
      for ( const auto& view : mFeature_infos[i].mView_list )
      {
        if ( true == mPic_query_mask[view.camera] )
        {
          auto mapIterImg = mapImgOrginToQuery.find( view.camera );
          if ( mapIterImg != mapImgOrginToQuery.end( ))
          {
            kCamQuery = mapIterImg->second;
            vPairImgFeatTo3DPointMatch[kCamQuery].push_back( std::make_pair( view.key, cntPointsInDB ) );
          }
          else{
            std::cout << "map find error. parsebundler.cpp line 646" << std::endl;
          }
        }
      }
      cntPointsInDB++;
    }
  }
  std::cout << "num of db points: " << cntPointsInDB << std::endl;

  // then save the match information
  std::ofstream os( sTrueMatchFile, std::ios::trunc );
  if ( !os.is_open() ){
    std::cout << "Open trueMatchFile fail: " << sTrueMatchFile << std::endl;
    return;
  }
  
  // rank the match
  size_t kMaximMatches = 0;
  for ( size_t kImgQuery = 0; kImgQuery < cntImgQuery; kImgQuery++ )
  {
    auto & vImgFeatMatch = vPairImgFeatTo3DPointMatch[kImgQuery];
    std::sort( vImgFeatMatch.begin(), vImgFeatMatch.end() );
    kMaximMatches = std::max( kMaximMatches, vImgFeatMatch.size() );
  }

  // save the match 
  size_t cntMatch= 0;
  for ( size_t kImgQuery = 0; kImgQuery < cntImgQuery; kImgQuery++ )
  {
    const auto & vImgFeatMatch = vPairImgFeatTo3DPointMatch[kImgQuery];
   
    // first line: img_id, feature_id1,...,feature_idn
    cntMatch = 0;
    os << kImgQuery << " ";
    for ( const auto& feat_match : vImgFeatMatch )
    { 
      os << feat_match.first << " ";
      ++cntMatch;
    }
    while ( cntMatch++ < kMaximMatches )
    { 
      os << -1 << " ";
    }
    os << std::endl;

    // second line: num_match, point_id1,...,point_idn
    cntMatch = 0;
    os << vImgFeatMatch.size() << " ";
    for ( const auto& feat_match : vImgFeatMatch )
    {
      os << feat_match.second << " ";
      ++cntMatch;
    }
    while ( cntMatch++ < kMaximMatches )
    {
      os << -1 << " ";
    }
    os << std::endl; 
  }
  os.close();
  //size_t kDBPoint = vPairImgFeatTo3DPointMatch[1][0].second;
  //for ( auto mapMember : map3DPointsOrigToDB ){
  //  if ( mapMember.second == kDBPoint ){
  //    std::cout << "original point index: " << mapMember.first << std::endl;
  //    for ( auto &viewPoint : mFeature_infos[mapMember.first].mView_list ){
  //      std::cout << viewPoint.camera << " " << viewPoint.key << std::endl;
  //    }
  //  }
  //}

  /* save the 3d point index map */
  // then save the match information
  std::ofstream ofs( sTrueMatchFile+".points.map", std::ios::trunc );
  if ( !ofs.is_open() ){
    std::cout << "Open points.map fail: " << sTrueMatchFile + ".points.map" << std::endl;
    return;
  }
  for (auto iter = map3DPointsOrigToDB.cbegin(); iter != map3DPointsOrigToDB.cend(); ++iter )
  { 
    ofs << iter->first << " " << iter->second << std::endl;
  }
  ofs.close();

}


void PARSE_BUNDLER::SetQueryMask(){
  mPic_query_mask.assign( mNumCameras, true );
}
void PARSE_BUNDLER::ResetQueryMask( ){
  mPic_query_mask.assign( mNumCameras, false );
}

//save the built information so next time directly load the file
//format::
//#cameras  #3d points bSaveRGB bSavedesc
//for each 3d points, contain the full information
//3dpoint(x, y, z, r g b)  
//#view  each view(camera, key x, y, scale, orientation, r, g, b)
//also all descriptor of one 3d points 
bool PARSE_BUNDLER::SaveFeature3DInfo( const std::string&s, bool bSaveRGB, bool bSavedesc) const
{
	if (mNumbPoints != mFeature_infos.size()){
		std::cout << "PARSE_BUNDLER: mNumbPoints != mFeature_infos.size() when save." << endl;
		return 0;
	}

	std::ofstream os(s, std::ios::trunc|std::ios::out);
	if (!os.is_open()){
		std::cout << "open parsed_bundler fail when save." << endl;
		return 0;
	}

  os << mNumCameras << " " << mNumbPoints << " " << (int)bSaveRGB << " " << (int)bSavedesc << std::endl;

	//each 3d points a line
	for (auto& feat_3d_info : mFeature_infos){
		os  << feat_3d_info.mPoint.x << " " 
        << feat_3d_info.mPoint.y << " " 
        << feat_3d_info.mPoint.z << " "
		    << (int)feat_3d_info.mPoint.r << " " 
        << (int)feat_3d_info.mPoint.g << " " 
        << (int)feat_3d_info.mPoint.b << " ";
		os  << std::endl;
		os  << feat_3d_info.mView_list.size() << " ";
		//save view
		for (auto& view : feat_3d_info.mView_list){
			os << view.camera << " " << view.key << " " << view.y << " "
				<< view.x << " " << view.scale << " " << view.orientation << " ";
      if ( bSaveRGB ){
        //std::cout << (int)view.rgb[0] << " " << (int)view.rgb[1] << " " << (int)view.rgb[2] << std::endl;
        os << (int)view.rgb[0] << " " << (int)view.rgb[1] << " " << (int)view.rgb[2] <<" ";
      }
		}
		os << std::endl;
    
    if ( !bSavedesc ) continue;
		//save decriptor num_desc = feat_3d_info.mDescriptor.size();
		for (auto& sift_desc : feat_3d_info.mDescriptor){
			for (int i = 0; i < sift_desc.legth; i++){
				os << int(sift_desc.ptrDesc[i]) << " ";
			}
		}
		os << std::endl;
	}

	os.close();
	return 1;
}

bool PARSE_BUNDLER::LoadFeature3DInfo(const std::string&s)
{
	std::ifstream is(s, std::ios::in);
	if (!is.is_open()){
		std::cout << " no parsed_bundler file then reload." << endl;
		return 0;
	}
  int bSaveRGB = 0, bSavedesc = 0;
  is >> mNumCameras >> mNumbPoints >> bSaveRGB >> bSavedesc;
	std::cout << "load feat_3d_info, total num of points:" << mNumbPoints << endl;
	
	mFeature_infos.clear();
	mFeature_infos.resize(mNumbPoints);
	
  int rgb[3];
  size_t view_list_len = 0;
	//each 3d points a line
	for (size_t i = 0; i < mNumbPoints; i++){
		auto& feat_3d_point = mFeature_infos[i].mPoint;
		is >> feat_3d_point.x >> feat_3d_point.y >> feat_3d_point.z;
    is >> rgb[0] >> rgb[1] >> rgb[2];
 
		view_list_len = 0;
		is >> view_list_len;

    feat_3d_point.r = unsigned char( rgb[0] );
    feat_3d_point.g = unsigned char( rgb[1] );
    feat_3d_point.b = unsigned char( rgb[2] );

		//load view_list
		mFeature_infos[i].mView_list.resize(view_list_len);
		for (auto& view : mFeature_infos[i].mView_list){
      // to be consistent with the bundler.out
			is >> view.camera >> view.key >> view.y >> view.x 
				>> view.scale >> view.orientation;
      if ( 1 == bSaveRGB ){
        is >> rgb[0] >> rgb[1] >> rgb[2];
        view.rgb[0] = unsigned char( rgb[0] );
        view.rgb[1] = unsigned char( rgb[1] );
        view.rgb[2] = unsigned char( rgb[2] );
      }
		}

    if ( 0 == bSavedesc ) continue;
		//load descriptor  num_desc = feat_3d_info.mDescriptor.size();
		mFeature_infos[i].mDescriptor.resize(view_list_len);
		for (int j = 0; j < view_list_len; j++){
			auto& sift_desc = mFeature_infos[i].mDescriptor[j];
			sift_desc.ptrDesc = new unsigned char[sift_desc.legth];
			if (!sift_desc.ptrDesc){
				std::cout << "new unsigned char fail. parsebundler.cpp line 330" << endl;
				j--;
				continue;
			}
			int tmp=0;
			for (int i = 0; i < sift_desc.legth; i++){
				is >> tmp; 
				sift_desc.ptrDesc[i] = unsigned char(tmp);
			}
		}
	}

	return 1;
}


std::vector< FEATURE_3D_INFO >& PARSE_BUNDLER::GetFeature3DInfo()
{
	return mFeature_infos;
}

const std::vector< FEATURE_3D_INFO >& PARSE_BUNDLER::GetFeature3DInfo() const
{
	return mFeature_infos;
}

ALL_PICTURES& PARSE_BUNDLER::GetAllPicturesAndCameras()
{
	return mAll_pic_cameras;
}

const ALL_PICTURES& PARSE_BUNDLER::GetAllPicturesAndCameras()const
{
	return mAll_pic_cameras;
}