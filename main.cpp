#include "config.h"
#include "parsebundler.h"

#ifdef USE_EXIF_READER
  #include "exif_reader/exif_reader.h"
#endif // USE_EXIF_READER

#include "opencv2/opencv.hpp"
#include <cstring>
#include <map>

enum eOperationMode
{
  BUNDLE_TO_PLY = 0,
  FIND_QUERY_TRUE_MATCH = 1,
  SCALE_3D_POINT_STAT = 2,
  BUNDLE_STATISTICS = 3,
  RESERVED = 4
};

int main(int argc, char** argv)
{
#if 0
  std::string imgfile1( "E:/Dubrovnik6K/query/33576593@N03_3126493232.jpg" );
  std::string imgfile2( "E:/Dubrovnik6K/query/croatia_kayak_2008_1888840167.jpg" );
  cv::Mat img1 = cv::imread( imgfile1, cv::IMREAD_COLOR );
  cv::Mat img2 = cv::imread( imgfile2, cv::IMREAD_COLOR );
  std::cout << " img1 rows: " << img1.rows << " cols: " << img1.cols 
    << " channel: " << img1.channels() << std::endl;
  std::cout << " img2 rows: " << img2.rows << " cols: " << img2.cols 
    << " channel: " << img2.channels() << std::endl;

  int img_width, img_height;
  exif_reader::open_exif( imgfile1.c_str() );
  img_width = exif_reader::get_image_width();
  img_height = exif_reader::get_image_height();
  std::cout << "img1 height: " << img_height << " width: " << img_width << std::endl;
  exif_reader::open_exif( imgfile1.c_str( ) );
  img_width = exif_reader::get_image_width();
  img_height = exif_reader::get_image_height();
  std::cout << "img2 height: " << img_height << " width: " << img_width << std::endl;

  std::vector<cv::Vec2d> feat_pts;
  feat_pts.push_back( cv::Vec2d( 487.31, 491.5 ) );
  feat_pts.push_back( cv::Vec2d( 580.37, 176.18 ) );
  feat_pts.push_back( cv::Vec2d( 676.69, 1037.89 ) );
  feat_pts.push_back( cv::Vec2d( 435.65, 412.84 ) );
  unsigned int row = 0 , col = 0;
  for ( auto point : feat_pts ){
    row = round( point[0] ), col = round( point[1] );
    std::cout << row << " " << col << " "
      << (int)img1.at<cv::Vec3b>( row, col )[2] << " "
      << (int)img1.at<cv::Vec3b>( row, col )[1] << " "
      << (int)img1.at<cv::Vec3b>( row, col )[0] << std::endl;

    row = floor( point[0] ), col = floor( point[1] );
    std::cout << row << " " << col << " "
      << (int)img1.at<cv::Vec3b>( row, col )[2] << " "
      << (int)img1.at<cv::Vec3b>( row, col )[1] << " "
      << (int)img1.at<cv::Vec3b>( row, col )[0] << std::endl;

    row = ceil( point[0] ), col = ceil( point[1] );
    std::cout << row << " " << col << " "
      << (int)img1.at<cv::Vec3b>( row, col )[2] << " "
      << (int)img1.at<cv::Vec3b>( row, col )[1] << " "
      << (int)img1.at<cv::Vec3b>( row, col )[0] << std::endl << std::endl;

  }

  //cv::imshow( "img1", img1 );
  //cv::imshow( "img2", img2 );
  //cv::waitKey();
  return 1;
#endif
#if 0
  // find image which should be adjusted.
  std::string imglist = "E:/Dubrovnik6K/list.db.txt";
  std::ifstream is( imglist, std::ios::_Nocreate );
  if ( !is.is_open() ){
    std::cout << "open img list fail: " << imglist << std::endl;
  }

  std::string line;
  std::vector<std::string> vImglist;
  vImglist.reserve( 800 );
  while ( getline( is, line ) ){
    std::istringstream newline( line );
    newline >> line;
    vImglist.push_back( line );
  }

  std::ofstream os( imglist + ".size" , std::ios::trunc);
  if ( !os.is_open() ){
    std::cout << "open image list.txt.size fail:" << imglist + ".size" << std::endl;
    return 0;
  }

  for ( size_t i = 0; i < vImglist.size(); ++i ){
    std::string imgfile( std::string( argv[1] ) +"/"+ vImglist[i] );
    cv::Mat img = cv::imread( imgfile );
    std::cout << "row: " << img.rows << " col: " << img.cols << std::endl;

    exif_reader::open_exif( imgfile.c_str( ) );
    int img_width, img_height;
    img_width = exif_reader::get_image_width();
    img_height = exif_reader::get_image_height();
    std::cout << "height: " << img_height << " width: " << img_width<< std::endl;
    if ( img.rows != img_height ){
      std::cout << "wrong size: " << imgfile << std::endl;
      os << vImglist[i] << std::endl;
    }
  }
  
  return 1;

#endif
  if ( argc < 5 ){
    std::cout << "___________________________________________________________" << std::endl;
    std::cout << "-                     bundler processor                   -" << std::endl;
    std::cout << "- Parameters:                                             -" << std::endl;
    std::cout << "- data_path  (have bundle in it)                          -" << std::endl;
    std::cout << "-     Bundler result data path.                           -" << std::endl;
    std::cout << "- bundle_out                                              -" << std::endl;
    std::cout << "-     bundler.out file name.                              -" << std::endl;
    std::cout << "- picture_list                                            -" << std::endl;
    std::cout << "-     picture_list_txt contain all pictures               -" << std::endl;
    std::cout << "- process                                                 -" << std::endl;
    std::cout << "-     what kind of process do with bundle.out             -" << std::endl;
    std::cout << "- 0: convert bundle.out to ply to be visualized           -" << std::endl;
    std::cout << "- 1: find the true match for query feature                -" << std::endl;
    std::cout << "- 2: do statics about scale of every 3d point             -" << std::endl;
    std::cout << "- 3: do statics about the color of every 3d point         -" << std::endl;
    std::cout << "- 4: ...                                                  -" << std::endl;
    std::cout << "-                                                         -" << std::endl;
    std::cout << "___________________________________________________________" << std::endl;
    return 0;
  }

  std::string data_path( argv[1] );
  std::string bundle_out( data_path + argv[2] );
  std::string ply_filename( bundle_out + ".ply" );
  std::string stat_result( bundle_out + ".stat" );
  std::string pic_list( argv[3] );

  std::cout << "bundle_out: " << bundle_out << std::endl;
  std::cout << "output file name: " << ply_filename << std::endl;
  std::cout << "output stat_res file name: " << stat_result << std::endl;

  eOperationMode eOperationmode;
  switch ( *argv[4] )
  {
    case '0': eOperationmode = BUNDLE_TO_PLY; break;
    case '1': eOperationmode = FIND_QUERY_TRUE_MATCH; break;
    case '2': eOperationmode = SCALE_3D_POINT_STAT; break;
    case '3': eOperationmode = BUNDLE_STATISTICS; break;
    default : eOperationmode = RESERVED; break;
  }

  PARSE_BUNDLER parse_bundller;
  parse_bundller.ParseBundlerFile( bundle_out );
 

  // convert the bundle file to ply so we can visualize the 3d point cloud with MeshLab 
  if ( eOperationmode == BUNDLE_TO_PLY ){
    std::cout << "operation mode: BUNDLE_TO_PLY" << std::endl;
    parse_bundller.ConvertPly( ply_filename );
  }
  // find the true match for the query image, so we can evaluate the match algorithm
  else if ( FIND_QUERY_TRUE_MATCH == eOperationmode ){
    std::cout << "operation mode: FIND_QUERY_TRUE_MATCH" << std::endl;
    std::string strListTxt( data_path + "/" + pic_list );
    parse_bundller.FindQueryPicture( strListTxt );
    parse_bundller.FindQueryFeatureTrueMatch( data_path + "/truematch.query" );

  }
  // color statistics
  else if ( eOperationmode == BUNDLE_STATISTICS ){
    std::cout << "operation mode: BUNDLE_STATISTICS" << std::endl;
    std::ifstream is( bundle_out + ".rgb", std::ios::_Nocreate );
    if ( is.is_open() ){
      std::cout << "LoadFeature3DInfo " << std::endl;
      parse_bundller.LoadFeature3DInfo( bundle_out + ".rgb" );
    }
    else{
      std::cout << "ReLoad Feature3DInfo " << std::endl;
      auto& pics = parse_bundller.GetAllPicturesAndCameras();
      pics.SetParameters( data_path, pic_list );
      pics.SetQueryFlag( true );
      pics.LoadPicturesKeyFile( false );
      parse_bundller.LoadFeatPointRGB();
      parse_bundller.LoadCameraInfo( false );
      std::cout << "load finish now save the feature_3d_info" << std::endl;
      parse_bundller.SaveFeature3DInfo( bundle_out + ".rgb", true, false );
    }

    parse_bundller.DoStatPoints( stat_result );
  }
  // scale statistics
  else if ( SCALE_3D_POINT_STAT == eOperationmode )
  {
    std::cout << "operation mode: SCALE_3D_POINT_STAT" << std::endl;
    auto& pics = parse_bundller.GetAllPicturesAndCameras();
    pics.SetParameters( data_path, pic_list );
    pics.SetQueryFlag( true );
    pics.LoadPicturesKeyFile( false );
    parse_bundller.LoadCameraInfo( false );

    std::ofstream os( bundle_out + ".scale", std::ios::trunc );
    size_t cntPoints = parse_bundller.GetNumPoints();
    const auto & feat_3d_infos = parse_bundller.GetFeature3DInfo();
    const auto & cameras = parse_bundller.GetAllPicturesAndCameras().GetAllCameras();
    
    for ( size_t cntPoint = 0; cntPoint < cntPoints; ++cntPoint )
    {
      cv::Vec3d point_in_world( feat_3d_infos[cntPoint].mPoint.x, feat_3d_infos[cntPoint].mPoint.y, feat_3d_infos[cntPoint].mPoint.z );
      
      os << point_in_world[0] << " " << point_in_world[1] << " " << point_in_world[2] << " "
         << feat_3d_infos[cntPoint].mView_list.size( ) << std::endl;
        for ( const auto & view : feat_3d_infos[cntPoint].mView_list )
        {
          cv::Vec3d point_in_frame = cameras[view.camera].rotation*point_in_world + cameras[view.camera].translation;
          cv::Vec3d cam_position = -cameras[view.camera].rotation.t( )*cameras[view.camera].translation;

          os << cam_position[0]   << " " << cam_position[1]   << " " << cam_position[2]   << " ";
          os << point_in_frame[0] << " " << point_in_frame[1] << " " << point_in_frame[2] << " ";
          os << view.camera << " " << view.key <<" "<< view.x << " " << view.y << " " << view.scale << " " << cameras[view.camera].focal_length << std::endl;
        }
      
    }
    os.close();
  }
  else std::cout << "operation mode: RESERVED" << std::endl;

  return 1;
}