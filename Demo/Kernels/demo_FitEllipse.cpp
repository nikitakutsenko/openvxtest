#include "../stdafx.h"

#include <opencv2/opencv.hpp>
#include <vector>

using namespace std;
using namespace cv;

extern "C"
{
#include <Kernels/ref.h>
#include <types.h>
}

#include "../DemoEngine.h"


class demo_FitEllipse : public IDemoCase
{
public:
   demo_FitEllipse(){}

   ///@see IDemoCase::ReplyName
   virtual std::string ReplyName() const override
   {
      return "FitEllipse";
   }

private:
   ///@see IDemoCase::execute
   virtual void execute() override;

   ///@brief provide interactive demo
   static void applyParameters(void* data);

private:
   cv::Mat m_srcImage;
};

///////////////////////////////////////////////////////////////////////////////
namespace
{
   const std::string m_openVXWindow    = "openVX";
   const std::string m_openCVWindow    = "openCV";
   const std::string m_originalWindow  = "original";
   const std::string m_diffWindow      = "Diff of " + m_openVXWindow + " and " + m_openCVWindow;
}

///////////////////////////////////////////////////////////////////////////////
void demo_FitEllipse::execute()
{

   cv::namedWindow(m_originalWindow,    CV_WINDOW_NORMAL);
   cv::namedWindow(m_openVXWindow,      CV_WINDOW_NORMAL);
   cv::namedWindow(m_openCVWindow,      CV_WINDOW_NORMAL);
   cv::namedWindow(m_diffWindow,        CV_WINDOW_NORMAL);

   const std::string imgPath = "../Image/apple.png";
   m_srcImage = cv::imread(imgPath, CV_LOAD_IMAGE_GRAYSCALE);
   cv::resize(m_srcImage, m_srcImage, cv::Size(500, 500), 0, 0, CV_INTER_LINEAR);
   cv::imshow(m_originalWindow, m_srcImage);

   applyParameters(this);
   cv::waitKey(0);
}

//////////////////////////////////////////////////////////////////////////////

void demo_FitEllipse::applyParameters(void* data) {
    auto demo = static_cast<demo_FitEllipse*>(data);
    const cv::Size imgSize(demo->m_srcImage.cols, demo->m_srcImage.rows);

    ///@{ OPENCV
    cv::Mat cvImage = cv::Mat::zeros(demo -> m_srcImage.size(), CV_8UC1), binImage;

    cv::Canny(demo -> m_srcImage, binImage, 50, 200, 3);
    vector<vector<Point>> contours;

    cv::findContours(binImage, contours, CV_RETR_TREE, CV_CHAIN_APPROX_SIMPLE);

	for (size_t i = 0; i < contours.size(); i++) {
        if (contours[i].size() >= 5) {
		    cv::RotatedRect box = fitEllipse(contours[i]);
            cv::ellipse(cvImage, box, Scalar(255, 255, 255));
        }
    }

    cv::imshow(m_openCVWindow, cvImage);

    ///@}


    ///@{ OPENVX
    cv::Mat vxImage = cv::Mat::zeros(demo->m_srcImage.size(), CV_8UC1);

    double scale = 1;

	for (size_t i = 0; i < contours.size(); i++) {
        if (contours[i].size() >= 5) {
            _vx_array box_, xs, ys;
            
            box_.size = 5;
            box_.data = (vx_float32*) calloc(5, sizeof(vx_float32));
            box_.array_type = VX_TYPE_FLOAT32;

            xs.size = contours[i].size();
            xs.array_type = VX_TYPE_FLOAT32;

            ys.size = contours[i].size();
            ys.array_type = VX_TYPE_FLOAT32;

            vx_float32 * _xs = new vx_float32[contours[i].size()];
            vx_float32 * _ys = new vx_float32[contours[i].size()];

            for (int j = 0; j < contours[i].size(); j++){
                _xs[j] = (vx_float32) contours[i][j].x / scale;
                _ys[j] = (vx_float32) contours[i][j].y / scale;
            }

            xs.data = _xs;
            ys.data = _ys;

		    vx_status st = ref_FitEllipse(&xs, &ys, &box_);
            if (st == VX_SUCCESS) {
                vx_float32 * bb = (vx_float32*) box_.data; 
                cv::RotatedRect box(
                    cv::Point2f(bb[0], bb[1]), 
                    cv::Size2f (bb[2], bb[3]), 
                    bb[4]
                );

                cv::ellipse(vxImage, box, Scalar(255, 255, 255));
            }
        }
    }

    cv::imshow(m_openVXWindow, vxImage);
    
    ///@}
    
    
    // Show difference of OpenVX and OpenCV
    const cv::Mat diffImage(imgSize, CV_8UC1);
    cv::absdiff(vxImage, cvImage, diffImage);
    cv::imshow(m_diffWindow, diffImage);
    //
}

///////////////////////////////////////////////////////////////////////////////
IDemoCasePtr CreateFitEllipseDemo() {
    return std::make_unique<demo_FitEllipse>();
}

