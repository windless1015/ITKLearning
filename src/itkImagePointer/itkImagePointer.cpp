
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkGDCMImageIO.h"
#include "itkImageToVTKImageFilter.h"
#include "itkGDCMImageIOFactory.h"
#include "itkNrrdImageIOFactory.h"
#include "itkNiftiImageIOFactory.h"
#include "itkMINCImageIOFactory.h"
#include "itkMetaImageIOFactory.h"
#include "itkBMPImageIOFactory.h"
#include "itkPNGImageIOFactory.h"
#include "itkJPEGImageIOFactory.h"
#include "itkGDCMSeriesFileNames.h"
#include "itkImageSeriesReader.h"
#include "itkOrientImageFilter.h"
#include "itkMetaDataObject.h"
#include "itkImageFileReader.h"

class VolumeImage
{
public:
    VolumeImage();
    ~VolumeImage();
    itk::Image<float, 3U>::Pointer getITKImagePointer() { return itkImage; }
    void readITKImage();
private:
    itk::Image<float, 3U>::Pointer itkImage;
};

VolumeImage::VolumeImage()
{
}

VolumeImage::~VolumeImage()
{
}

void VolumeImage::readITKImage()
{

}


int main(int argc, char * argv[])
{
    VolumeImage itkVolumeImage;

    /*
    using PixelType = float;
    using ImageType = itk::Image<PixelType, 3>;
    //1. read the image from local file
	itk::ImageFileReader<ImageType>::Pointer reader = itk::ImageFileReader<ImageType>::New();
	reader->SetFileName(path);
	try {
		reader->Update();
	}
	catch (itk::ExceptionObject & err) {
		std::cout << "ERROR while parsing fixed volume: " << err.what() << std::endl;
		return false;
	}

	//2. OrientImageFilter change image orientation(direction)
	itk::OrientImageFilter<ImageType, ImageType>::Pointer orientationFilter = itk::OrientImageFilter<ImageType, ImageType>::New();
	orientationFilter->UseImageDirectionOn();
	orientationFilter->SetDesiredCoordinateOrientation(itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_RAI);
	orientationFilter->SetInput(reader->GetOutput());
	try {
		orientationFilter->Update();
	}
	catch (itk::ExceptionObject & err) {
		std::cout << "ERROR while correcting fixed volume orientation: " << err.what() << std::endl;
		return false;
	}

	//3. get image that has been changed orientation
	itk::Image<float, 3U> *volumeImgPtr = orientationFilter->GetOutput();
    */

    //std::string filePath = "D:\\projects\\Fussen\\shangtaoSkull\\dicom\\0262.dcm";

    //using PixelType = signed short;
    //using ImageType = itk::Image<PixelType, 3>;
    ////1. read the image from local file
    //itk::ImageFileReader<ImageType>::Pointer reader = itk::ImageFileReader<ImageType>::New();
    //reader->SetFileName(filePath);
    //try {
    //    reader->Update();
    //}
    //catch (itk::ExceptionObject & err) {
    //    std::cout << "ERROR " << err.what() << std::endl;
    //    return false;
    //}

    //itk::OrientImageFilter<ImageType, ImageType>::Pointer orientationFilter = itk::OrientImageFilter<ImageType, ImageType>::New();
    //orientationFilter->UseImageDirectionOn();
    //orientationFilter->SetDesiredCoordinateOrientation(itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_RAI);
    //orientationFilter->SetInput(reader->GetOutput());
    //try {
    //    orientationFilter->Update();
    //}
    //catch (itk::ExceptionObject & err) {
    //    std::cout << "ERROR while correcting fixed volume orientation: " << err.what() << std::endl;
    //    return false;
    //}

    ////3. get image that has been changed orientation
    ////itk::Image<PixelType, 3U> *volumeImgPtr = orientationFilter->GetOutput();

    ////创建读取DCM的GDCMIOImage类
    //using ImageIOType = itk::GDCMImageIO;
    //ImageIOType::Pointer gdcmImageIO = ImageIOType::New();
    //reader->SetImageIO(gdcmImageIO);
    //try
    //{
    //    reader->Update();
    //}
    //catch (itk::ExceptionObject& e)
    //{
    //    std::cerr << "exception in file reader" << std::endl;
    //    std::cerr << e << std::endl;
    //    return EXIT_FAILURE;
    //}

    ////链接滤波器，转换为VTK类型
    //typedef itk::ImageToVTKImageFilter<ImageType> ConnectorType;
    //ConnectorType::Pointer connector = ConnectorType::New();
    //connector->SetInput(reader->GetOutput());
    //try
    //{
    //    connector->Update();
    //}
    //catch (itk::ExceptionObject& e)
    //{
    //    std::cerr << "exception in file reader" << std::endl;
    //    std::cerr << e << std::endl;
    //    return EXIT_FAILURE;
    //}

    ///* vtkImageActor在3D场景下渲染图像 */
    //vtkSmartPointer<vtkImageActor> actor = vtkImageActor::New();
    //actor->SetInputData(connector->GetOutput());
    //actor->InterpolateOff();
    //actor->Update();


    //vtkSmartPointer<vtkRenderer> render = vtkRenderer::New();
    //render->AddActor(actor);
    //render->SetBackground(255, 1, 1);//设置窗口背景颜色

    //vtkSmartPointer<vtkRenderWindow> window = vtkRenderWindow::New();
    //window->SetSize(600, 600);       //设置窗口大小
    //window->AddRenderer(render);
    //window->SetWindowName("the CT picture ");//设置窗口名称

    //vtkSmartPointer<vtkRenderWindowInteractor> interactor = vtkRenderWindowInteractor::New();
    //interactor->SetRenderWindow(window);

    //interactor->Initialize();
    //interactor->Start();
    //return 0;

    return EXIT_SUCCESS;


}

