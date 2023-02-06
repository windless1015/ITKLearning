//#include "itkImage.h"
//#include "itkGDCMImageIO.h"
//#include "itkGDCMSeriesFileNames.h"
//#include "itkImageSeriesReader.h"
//#include "itkImageFileWriter.h"
#include <vector>
#include <string>
#include <algorithm>
#include <iostream>

#include <vtkSmartPointer.h>
#include <vtkNrrdReader.h>

#include <itkImage.h>
#include <itkImageSeriesReader.h>
#include <itkGDCMImageIO.h>
#include <itkGDCMSeriesFileNames.h>
#include <itkMacro.h>
#include <itkImageRegionIteratorWithIndex.h>
#include <itkImageFileWriter.h>
#include <itkImageFileReader.h>
#include <itkMetaDataObject.h>
#include <itkNrrdImageIO.h>
#include <itkVectorImage.h>
#include <itkVariableLengthVector.h>
#include <itkOrientImageFilter.h>

#include "itkGDCMImageIOFactory.h"
#include "itkNrrdImageIOFactory.h"
#include "itkNiftiImageIOFactory.h"
#include "itkMINCImageIOFactory.h"
#include "itkMetaImageIOFactory.h"
#include "itkVTKImageToImageFilter.h"

//Write binary files to disk, with Number of elements.

void writeToBin(float *Output, int Num_Elements, const std::string FILENAME) {
    std::ofstream OutputStream;
    OutputStream.open(FILENAME, std::ios::app | std::ios::binary);

    if (!OutputStream.good()) {
        //std::cerr("Failed to open " + FILENAME);
        exit(0);
    }

    OutputStream.write(reinterpret_cast<char*>(Output), sizeof(float) * Num_Elements);

    OutputStream.close();
}

void readInformation(std::string filePath, int* dimensions, float* origins, double* directions, double* spacing, std::vector<float>& buffer);

bool writeNrrdFromITKWiki(std::string, std::string);
void writeBySelf(std::string outputPath, int* dimensions, float* origins, double* directions, double* spacing, std::vector<float>& buffer);
void readNrrdImageAndGenerateNrrd(std::string srcFile, std::string outFile);

int main(int argc, char * argv[])
{
    //mha
    std::string path1 = "D:\\Ultrast\\Patients\\RD_UT20220316143157\\appt_1\\volume_datasets\\3D-USScan_1.mha";
    std::string path2 = "D:\\Ultrast\\Patients\\RD_UT20220316143157\\appt_1\\volume_datasets\\MainMR_t1_tse_tra.mha";
    std::string pathNrrd = "D:\\Ultrast\\Patients\\RD_UT20220316143157\\appt_1\\volume_datasets\\fusion\\3D-Scan_1.nrrd";
    //this function is testing the case of converting vtkImageData into nrrd type
    readNrrdImageAndGenerateNrrd("D:/test.nrrd", "D:/tttttt.nrrd");

    //this function is testing the case of generating the nrrd file by myself
    int dimensions[3];
    float origins[3];
    double directions[3];
    double spacing[3];

    std::vector<float> buffer;
    readInformation(path1, dimensions, origins, directions, spacing, buffer);
    writeBySelf("D:/generateNrrdByself.nrrd", dimensions, origins, directions, spacing, buffer);



    //this function is testing the case of generating the nrrd file by official method
    writeNrrdFromITKWiki(path1, "D:/officialMethodNrrd.nrrd");
    return EXIT_SUCCESS;

}

//this function comes from https://examples.itk.org/src/io/gdcm/readdicomseriesandwrite3dimage/documentation
bool writeNRRDFromITKSamples()
{
    //std::string dirName = "."; // current directory by default

    //using PixelType = float;
    //constexpr unsigned int Dimension = 3;
    //using ImageType = itk::Image<PixelType, Dimension>;

    //using NamesGeneratorType = itk::GDCMSeriesFileNames;
    //auto nameGenerator = NamesGeneratorType::New();

    //nameGenerator->SetUseSeriesDetails(true);
    //nameGenerator->AddSeriesRestriction("0008|0021");
    //nameGenerator->SetGlobalWarningDisplay(false);
    //nameGenerator->SetDirectory(dirName);

    //try
    //{
    //    using SeriesIdContainer = std::vector<std::string>;
    //    const SeriesIdContainer & seriesUID = nameGenerator->GetSeriesUIDs();
    //    auto                      seriesItr = seriesUID.begin();
    //    auto                      seriesEnd = seriesUID.end();

    //    if (seriesItr != seriesEnd)
    //    {
    //        std::cout << "The directory: ";
    //        std::cout << dirName << std::endl;
    //        std::cout << "Contains the following DICOM Series: ";
    //        std::cout << std::endl;
    //    }
    //    else
    //    {
    //        std::cout << "No DICOMs in: " << dirName << std::endl;
    //        return EXIT_SUCCESS;
    //    }

    //    while (seriesItr != seriesEnd)
    //    {
    //        std::cout << seriesItr->c_str() << std::endl;
    //        ++seriesItr;
    //    }

    //    seriesItr = seriesUID.begin();
    //    while (seriesItr != seriesUID.end())
    //    {
    //        std::string seriesIdentifier;
    //        if (argc > 3) // If seriesIdentifier given convert only that
    //        {
    //            seriesIdentifier = argv[3];
    //            seriesItr = seriesUID.end();
    //        }
    //        else // otherwise convert everything
    //        {
    //            seriesIdentifier = seriesItr->c_str();
    //            seriesItr++;
    //        }
    //        std::cout << "\nReading: ";
    //        std::cout << seriesIdentifier << std::endl;
    //        using FileNamesContainer = std::vector<std::string>;
    //        FileNamesContainer fileNames = nameGenerator->GetFileNames(seriesIdentifier);

    //        using ReaderType = itk::ImageSeriesReader<ImageType>;
    //        auto reader = ReaderType::New();
    //        using ImageIOType = itk::GDCMImageIO;
    //        auto dicomIO = ImageIOType::New();
    //        reader->SetImageIO(dicomIO);
    //        reader->SetFileNames(fileNames);
    //        reader->ForceOrthogonalDirectionOff(); // properly read CTs with gantry tilt

    //        std::string outFileName;
    //        if (argc > 2)
    //        {
    //            outFileName = argv[2];
    //        }
    //        else
    //        {
    //            outFileName = dirName + std::string("/") + seriesIdentifier + ".nrrd";
    //        }
    //        std::cout << "Writing: " << outFileName << std::endl;
    //        try
    //        {

    //            itk::WriteImage(reader->GetOutput(), outFileName, true); // compression
    //        }
    //        catch (const itk::ExceptionObject & ex)
    //        {
    //            std::cout << ex << std::endl;
    //            continue;
    //        }
    //    }
    //}
    //catch (const itk::ExceptionObject & ex)
    //{
    //    std::cout << ex << std::endl;
    //    return EXIT_FAILURE;
    //}
    return EXIT_SUCCESS;
}

bool writeNrrdFromITKWiki(std::string srcFile, std::string outFile)
{
    //typedef signed short                      PixelType;
    typedef float                      PixelType;
    typedef itk::VectorImage<PixelType, 3>	    DiffusionImageType;
    typedef DiffusionImageType::Pointer				DiffusionImagePointer;


    typedef itk::ImageFileReader<DiffusionImageType,
        itk::DefaultConvertPixelTraits< PixelType > > FileReaderType;
    FileReaderType::Pointer reader = FileReaderType::New();
    reader->SetFileName(srcFile);
    reader->Update();

    itk::NrrdImageIO::Pointer io = itk::NrrdImageIO::New();
    //io->SetNrrdVectorType(nrrdKindList);
    //io->set
    io->SetFileType(itk::ImageIOBase::Binary);

    typedef itk::ImageFileWriter<DiffusionImageType > WriterType;
    WriterType::Pointer nrrdWriter = WriterType::New();
    nrrdWriter->UseInputMetaDataDictionaryOn();
    nrrdWriter->SetInput(reader->GetOutput());
    nrrdWriter->SetImageIO(io);
    nrrdWriter->SetFileName(outFile);
    try
    {
        nrrdWriter->Update();
    }
    catch (itk::ExceptionObject e)
    {
        std::cout << e << std::endl;
    }

    return 0;
}



void writeBySelf(std::string outputPath, int* dimensions, float* origins, double* directions, double* spacing, std::vector<float>& buffer)
{
    std::ofstream Export_NRRD;
    Export_NRRD.open(outputPath, std::ofstream::out);

    Export_NRRD << "NRRD0004" << std::endl;
    Export_NRRD << "type: float" << std::endl;  //"type: uint8"

    Export_NRRD << "dimension: 3" << std::endl;
    /*Export_NRRD << "space: scanner-xyz" << std::endl;*/  //refer to https://teem.sourceforge.net/nrrd/format.html 
    Export_NRRD << "space: left-posterior-superior" << std::endl;
    Export_NRRD << "sizes: " << dimensions[0] << " " << dimensions[1] << " " << dimensions[2] << std::endl;
    Export_NRRD << "space directions: (" << directions[0] << ",0,0) (0," << directions[1] << ",0) (0,0," << directions[2] << ")" << std::endl;
    Export_NRRD << "kinds: domain domain domain" << std::endl;
    Export_NRRD << "endian: little" << std::endl;
    Export_NRRD << "encoding: raw" << std::endl;
    Export_NRRD << "space origin: (" << origins[0] << "," << origins[1] << "," << origins[2] << ")" << std::endl << std::endl;
    Export_NRRD.close();

    writeToBin(buffer.data(), buffer.size(), outputPath);
}

void readInformation(std::string filePath, int* dimensions, float* origins, double* directions, double* spacing, std::vector<float>& dataBuffer)
{
    //extract the extension to determine what data type it is
    itk::NrrdImageIOFactory::RegisterOneFactory();
    itk::NiftiImageIOFactory::RegisterOneFactory();
    itk::MINCImageIOFactory::RegisterOneFactory();
    itk::MetaImageIOFactory::RegisterOneFactory();
    itk::GDCMImageIOFactory::RegisterOneFactory();

    using PixelType = float;
    using ImageType = itk::Image<PixelType, 3>;
    //1. read the image from local file
    itk::ImageFileReader<ImageType>::Pointer reader = itk::ImageFileReader<ImageType>::New();
    reader->SetFileName(filePath);
    try {
        reader->Update();
    }
    catch (itk::ExceptionObject & e) {
        std::cerr << "Volume Image -- ITK image filter error" << std::endl;
        return;
    }
    //itk::Image<ImageType>::Pointer itkImagePointer = reader->GetImageIO();
    itk::ImageIOBase* io = reader->GetImageIO();
    for (int i = 0; i < 3; i++)
    {
        dimensions[i] = io->GetDimensions(i);
        origins[i] = io->GetOrigin(i);
        directions[i] = io->GetDirection(i)[i];
        spacing[i] = io->GetSpacing(i);
    }
    
    //2. filer the orientation 
    itk::OrientImageFilter<ImageType, ImageType>::Pointer orientationFilter = itk::OrientImageFilter<ImageType, ImageType>::New();
    orientationFilter->UseImageDirectionOn();
    // RAI means Right to left, A means anterior to posterior, I means interior to superior, the same orientation in ITK called RAI, in DICOM called LPS, in Nrrd called left-posterior-superior
    orientationFilter->SetDesiredCoordinateOrientation(itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_RAI);
    orientationFilter->SetInput(reader->GetOutput());
    try {
        orientationFilter->Update();

        itk::Image<float, 3U>::Pointer volumeImgPtr = orientationFilter->GetOutput();
        int bufferSize = volumeImgPtr->GetPixelContainer()->Size();
        dataBuffer.resize(bufferSize);

        std::copy(volumeImgPtr->GetPixelContainer()->GetBufferPointer(),
            volumeImgPtr->GetPixelContainer()->GetBufferPointer() + bufferSize,
            dataBuffer.data());
    }
    catch (itk::ExceptionObject & e) {
    }

}


//read nrrd volume image with vtkNrrdReader and generate nrrd again
// in some cases, we have the vtkImageData already and we need to generate it into nrrd type
void readNrrdImageAndGenerateNrrd(std::string srcFile, std::string outFile)
{
    vtkSmartPointer<vtkNrrdReader> reader1 = vtkSmartPointer<vtkNrrdReader>::New();
    if (!reader1->CanReadFile(srcFile.c_str()))
    {
        std::cerr << "Reader reports " << srcFile << " cannot be read.";
        return;
    }
    reader1->SetFileName(srcFile.c_str());
    reader1->Update();

    using PixelType = float;
    using ImageType = itk::Image<PixelType, 3>;
    using FilterType = itk::VTKImageToImageFilter<ImageType>;
    auto filter = FilterType::New();
    filter->SetInput(reader1->GetOutput());

    try
    {
        filter->Update();
    }
    catch (const itk::ExceptionObject & error)
    {
        std::cerr << "Error: " << error << std::endl;
        return;
    }
    ImageType::ConstPointer itkImage = filter->GetOutput();



    itk::NrrdImageIO::Pointer io = itk::NrrdImageIO::New();
    io->SetFileType(itk::ImageIOBase::Binary);

    typedef itk::ImageFileWriter<ImageType > WriterType;
    WriterType::Pointer nrrdWriter = WriterType::New();
    nrrdWriter->UseInputMetaDataDictionaryOn();
    nrrdWriter->SetInput(itkImage);
    nrrdWriter->SetImageIO(io);
    nrrdWriter->SetFileName(outFile);
    try
    {
        nrrdWriter->Update();
    }
    catch (itk::ExceptionObject e)
    {
        std::cout << e << std::endl;
    }

}