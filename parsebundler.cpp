
#include "parsebundler.h"

#include "opencv2/opencv.hpp"

#include <string.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <utility>
#include <numeric>

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

  /* Print the ply header */
  fprintf( f, ply_header, mNumbPoints );
  for ( const auto& point_info : mFeature_infos ){
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
    float rgb_mean[3];
    float rgb_variance[3];
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
	string line;
	while (is >> line)
	{
		//stringstream istrstream;
		mPic_query_mask.push_back( line[0] == 'q' );
	}

}

//
void PARSE_BUNDLER::WriteQueryBundler( const std::string& s, bool bWritepoint ) const
{
  ofstream os( s, std::ios::out | std::ios::trunc );
  if ( 0 == os.is_open() ){
    std::cout << " open query bundle file fail: " << s << endl;
    return;
  }

  const auto & mPictures = mAll_pic_cameras.GetAllPictures();
  const auto & mCameras = mAll_pic_cameras.GetAllCameras();

  os << "# Bundle file v0.3" << endl;

  size_t num_cam = 0, num_pts = 0;
  //count to number of cameras
  for ( size_t i = 0; i < mPic_query_mask.size(); i++ ){
    if ( 1 == mPic_query_mask[i] )
    {
      ++num_cam;
    }
  }
  // find 3d points see in query pictures
  std::vector<bool> points_in_query( mNumbPoints, false );
  for ( size_t i = 0; bWritepoint && i < mNumbPoints; i++ )
  {
    for ( const auto& view : mFeature_infos[i].mView_list )
    {
      if ( 1 == mPic_query_mask[view.camera] ){
        points_in_query[i] = true;
        ++num_pts;
        break;
      }
    }
  }
  os << num_cam << " " << num_pts << endl;

  for ( size_t i = 0; i < mPic_query_mask.size(); i++ ){
    if ( 1 == mPic_query_mask[i] ){
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


  for ( size_t i = 0; bWritepoint && i < mNumbPoints; i++ )
  {
    if ( 1 == points_in_query[i] ){
      int view_lenth = mFeature_infos[i].mView_list.size();
      int view_len_query = 0;

      //calculate true view_length include only query points in cameras
      for ( int j = 0; j < view_lenth; j++ ){
        if ( 1 == mPic_query_mask[mFeature_infos[i].mView_list[j].camera] )
        {
          view_len_query++;
        }
      }

      os << mFeature_infos[i].mPoint.x << " "
        << mFeature_infos[i].mPoint.y << " "
        << mFeature_infos[i].mPoint.z << std::endl
        << (int)mFeature_infos[i].mPoint.r << " "
        << (int)mFeature_infos[i].mPoint.g << " "
        << (int)mFeature_infos[i].mPoint.b << std::endl
        << view_len_query << " ";

      for ( int j = 0; j < view_lenth; j++ )
      {
        if ( 1 == mPic_query_mask[mFeature_infos[i].mView_list[j].camera] )
        {
          os << mFeature_infos[i].mView_list[j].camera << " "
            << mFeature_infos[i].mView_list[j].key << " "
            << mFeature_infos[i].mView_list[j].x << " "
            << mFeature_infos[i].mView_list[j].y << " ";
        }
      }
      os << std::endl;
    }
  }

  os.close();
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
			os << view.camera << " " << view.key << " " << view.x << " "
				<< view.y << " " << view.scale << " " << view.orientation << " ";
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
	std::ifstream is(s, std::ios::_Nocreate);
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
			is >> view.camera >> view.key >> view.x >> view.y 
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
