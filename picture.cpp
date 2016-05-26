#include "picture.h"
#include "timer/timer.h"
#include "exif_reader/exif_reader.h"

#include <fstream>
#include <sstream>
#include <memory>
#include <omp.h>

using std::cout;
using std::endl;
using std::string;


int SIFT_Descriptor::legth = 128;


//calculate the distance between the two sift descriptor
int CalculateSIFTDistanceSquared(const unsigned char* d1, const unsigned char* d2)
{
	int dif, distsq = 0;
	for (int i = 0; i < 128; i++){
		dif = d1[i] - d2[i];
		distsq += dif*dif;
	}
	return distsq;
}

/****** for class PICTURE ****/
//default constructor
PICTURE::PICTURE(): mDes_length(0), mKeypoint_num(0), mImageHeight(0), mImgaeWidth(0){
	
}

PICTURE::PICTURE( const int height, const int width ) : 
  mDes_length( 0 ), mKeypoint_num( 0 ), 
  mImageHeight( height ), mImgaeWidth( width )
{
  
}

PICTURE::~PICTURE(){

};

//set and get image size
void PICTURE::SetImageSize(const int height, const int width)
{
	mImgaeWidth = width;
	mImageHeight = height;
}

//get size
void PICTURE::GetImageSize(int& height, int& width) const
{
	height = mImageHeight;
	width = mImgaeWidth;
}

//clear all the data, delete the pointer content
void PICTURE::ClearData()
{
	mDescriptors.clear();
	mFeature_points.clear();
}

//load  key points and descriptors  
bool PICTURE::LoadKeyPointAndDes(const std::string& des_filename, bool bLoadDesc)
{
	std::ifstream infile(des_filename, std::ios::in);
	if (!infile.is_open()){
		std::cout << " key file open fail" <<des_filename <<endl;
		return 0;
	}
  //std::cout << "LoadKeyPointAndDes: " << des_filename << " " << bLoadDesc << endl;
    
	infile >> mKeypoint_num >> mDes_length;
	assert(mKeypoint_num >= 0 && mDes_length >= 0);

	mFeature_points.resize(mKeypoint_num);
	mDescriptors.resize(mKeypoint_num);

  std::string line;
	//assert(mDes_length == mDescriptors[0].legth);
	for (int cnt = 0; cnt < mKeypoint_num; cnt++)
	{
		auto& sift_keypt = mFeature_points[cnt];
		// y x  scale  orientation
		infile >> sift_keypt.y >> sift_keypt.x 
			>> sift_keypt.scale >> sift_keypt.orientation;

    // now do not load desc
    if ( !bLoadDesc )
    {
      for ( int i = 0; i < 8; i++ )
      { 
        getline( infile, line ); 
      }
      continue;
    }

		//directly operate on the desc
		//do not use temp variable and then push_back it into the vector
		//since program end the temp variable will release the newed memory 
		auto& sift_desc = mDescriptors[cnt];

		sift_desc.ptrDesc = new unsigned char[sift_desc.legth];
		if (sift_desc.ptrDesc == nullptr || mDes_length != sift_desc.legth){
			std::cout << "new error(picture.cpp, line 06)" << std::endl;
			return 0;
		}

		//read the descriptor
		int des_temp = 0;
		for (int i = 0; i < sift_desc.legth; ++i)
		{
			infile >> des_temp;
			sift_desc.ptrDesc[i] = (unsigned char)des_temp;
		}

	}

	infile.close();
	return 1;
}

/*************************************/
// after load key point, also read the rgb of the feature point
bool PICTURE::LoadKeyPointRGB( const std::string& des_filename )
{
  //std::cout << "read rgb of feature point: " << std::endl;
  std::string img_file( des_filename );
  img_file.replace( img_file.end() - 3, img_file.end(), "jpg" );

  //std::cout << "LoadKeyPointRGB : " << img_file << std::endl;
  cv::Mat img = cv::imread( img_file, cv::IMREAD_COLOR );
  if ( img.empty()){
    std::cout << "open img fail: " << img_file << std::endl;
    return 0;
  }

  for ( auto& feat_pt : mFeature_points ){
    int row = round(feat_pt.y), col = round(feat_pt.x);
    if ( img.rows == mImageHeight && img.cols == mImgaeWidth 
      && row < mImageHeight && col < mImgaeWidth){
      feat_pt.rgb[0] = img.at<cv::Vec3b>( row, col )[2];
      feat_pt.rgb[1] = img.at<cv::Vec3b>( row, col )[1];
      feat_pt.rgb[2] = img.at<cv::Vec3b>( row, col )[0];
      //std::cout << int( feat_pt.rgb[0] ) << " " << int( feat_pt.rgb[1] ) << " " << int( feat_pt.rgb[2] ) << std::endl;
    }
    else{
      std::cout << "img size error:  height  width  row  col  row_pt  col_pt:" << std::endl;
      std::cout << mImageHeight << " " << mImgaeWidth << " " 
        << img.rows << " " << img.cols << " " << row << " " << col << std::endl;
    }  
  }

  return 1;
}



/*************************************/

/****** for Class ALL_PICTURES ***/
//set all empty
ALL_PICTURES::ALL_PICTURES() :mKeyfilepath(""), 
mImagepath(""), mPicturelistfile(""), mIsqueryimage(false){
	
}

//set key and image have the same path
ALL_PICTURES::ALL_PICTURES(const std::string& key, const std::string &list)
: mKeyfilepath(key), mImagepath(key), mPicturelistfile(list), mIsqueryimage(false){

}

//set separated path 
ALL_PICTURES::ALL_PICTURES(const std::string& key, const std::string& image, const std::string &list)
: mKeyfilepath(key), mImagepath(image), mPicturelistfile(list), mIsqueryimage(false){

}

//destructor
ALL_PICTURES::~ALL_PICTURES(){

}

//if there already exists pictures, clear them then reload
bool ALL_PICTURES::LoadPicturesKeyFile(bool bLoaddesc)
{
	Timer timer;
	timer.Start();

	//first clear all pictures if already loaded
	ClearPics();
  
	std::ifstream infile(mKeyfilepath + '/' + mPicturelistfile, std::ios::in);
	if (!infile.is_open()){
		std::cout << "Open list file fail: " << mPicturelistfile << endl;
		return 0;
	}

  //clear the pictures and then load pictures;
  mPictures.clear();
  mPictures.reserve( 1000 );
  mQueryfocal.reserve( 1000 );

	std::vector<string> pic_keyfilename;
	std::string  line_in_file;
	std::string	picture_filename;
  int	height = 0, width = 0, temp_zero = 0;
  double	temp_f = 0.0;

  while ( getline( infile, line_in_file ) )
  {
    std::istringstream words_in_line( line_in_file );
    words_in_line >> picture_filename >> height >> width;

    pic_keyfilename.push_back( picture_filename );
    mPictures.push_back( PICTURE(height, width) );
    //std::cout << "img " << picture_filename << " size: " << height << " "<<width << endl;
    if ( false == mIsqueryimage ) continue;

    //check if have f 
    if ( ( words_in_line >> temp_zero ) && temp_zero == 0 ){
      words_in_line >> temp_f;
      mQueryfocal.push_back( temp_f );
    }
    else mQueryfocal.push_back( -1.0 );
  }

  // now add the image size in the img list
#if 0
	//for query images, also load the image size;
	if (mIsqueryimage){
		// first we need to get the dimensions of the image
		// then center the keypoints around the center of the image
		std::string img_filename;
		for (size_t i = 0; i < mPictures.size(); i++)
		{
			img_filename = mKeyfilepath + '/' + pic_keyfilename[i];
			exif_reader::open_exif(img_filename.c_str());
			int img_width, img_height;
			img_width = exif_reader::get_image_width();
			img_height = exif_reader::get_image_height();
      cout << "image height: " << img_height << " ,width " << img_width << std::endl;
			mPictures[i].SetImageSize(exif_reader::get_image_height(), exif_reader::get_image_width());
			exif_reader::close_exif();
		}
	}
#endif

  std::cout << "load img file name, then load key" << std::endl;

#pragma omp parallel for
	for (int i = 0; i < pic_keyfilename.size(); i++)
	{
		pic_keyfilename[i].replace(pic_keyfilename[i].end() - 3, pic_keyfilename[i].end(), "key");
		//load a picture
    mPictures[i].LoadKeyPointAndDes( mKeyfilepath + "/" + pic_keyfilename[i], bLoaddesc );
    mPictures[i].LoadKeyPointRGB( mKeyfilepath + "/" + pic_keyfilename[i] );
  }

	timer.Stop();
	std::cout << "load picture keys time: " << timer.GetElapsedTimeAsString() << endl;
	return 1;
}


//if it is for query image, make flag true
void ALL_PICTURES::SetQueryFlag(const bool flag)
{
	mIsqueryimage = flag;
}

const bool ALL_PICTURES::RetQueryFlag() const
{
	return mIsqueryimage;
}

//set the string contents
void ALL_PICTURES::SetParameters(const std::string& model_path, const std::string &list)
{
	mKeyfilepath = model_path;
	mPicturelistfile = list;
}

//load the camera pose ground truth from bundler.query.out
bool ALL_PICTURES::LoadCamerasPose(const std::string& s)
{
	Timer timer;
	timer.Start();

	std::ifstream instream(s, std::ios::in);
	if (!instream.is_open()){
		std::cout << "open bundler fail: " << s << std::endl;
		return 0;
	}
	std::string line;
	std::getline(instream, line);//header

	int num_cam = 0, num_points = 0;
	instream >> num_cam >> num_points;
	mCameras.resize(num_cam);

	assert(mCameras.size()==mPictures.size());

	for (size_t i = 0; i < num_cam; i++){
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
		
		mPictures[i].GetImageSize(mCameras[i].height, mCameras[i].width);
	}

	instream.close();

	timer.Stop();
	std::cout << "load camera true pose time: " << timer.GetElapsedTimeAsString() << endl;
	return 1;
}