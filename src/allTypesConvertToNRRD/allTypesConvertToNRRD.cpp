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
#include "itkMetaDataDictionary.h"
#include <vtkMetaImageReader.h>

typedef short InputPixelType;
typedef float OutputPixelType;

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
void addMetaTagsToNRRDFile(std::string srcFile, std::string outFile);
int main(int argc, char * argv[])
{
    //mha
    std::string path1 = "D:\\Ultrast\\Patients\\RD_UT20220316143157\\appt_1\\volume_datasets\\3D-USScan_1.mha";
    std::string path2 = "D:\\Ultrast\\Patients\\RD_UT20220316143157\\appt_1\\volume_datasets\\MainMR_t1_tse_tra.mha";
    std::string pathNrrd = "D:\\Ultrast\\Patients\\RD_UT20220316143157\\appt_1\\volume_datasets\\fusion\\3D-Scan_1.nrrd";
    //this function is testing the case of converting vtkImageData into nrrd type
    //readNrrdImageAndGenerateNrrd(path1, "D:/tttttt.nrrd");
    addMetaTagsToNRRDFile("D:/tttttt.nrrd", "D:/tttmodified.nrrd");

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
    typedef itk::VectorImage<InputPixelType, 3>	    DiffusionImageType;
    typedef DiffusionImageType::Pointer				DiffusionImagePointer;


    typedef itk::ImageFileReader<DiffusionImageType,
        itk::DefaultConvertPixelTraits<InputPixelType>> FileReaderType;
    FileReaderType::Pointer reader = FileReaderType::New();
    reader->SetFileName(srcFile);
    reader->Update();

    itk::NrrdImageIO::Pointer io = itk::NrrdImageIO::New();
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

    using ImageType = itk::Image<InputPixelType, 3>;
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

        itk::Image<InputPixelType, 3U>::Pointer volumeImgPtr = orientationFilter->GetOutput();
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
    /*vtkSmartPointer<vtkNrrdReader> nrrdReader = vtkSmartPointer<vtkNrrdReader>::New();
    if (!nrrdReader->CanReadFile(srcFile.c_str()))
    {
        std::cerr << "Reader reports " << srcFile << " cannot be read.";
        return;
    }
    nrrdReader->SetFileName(srcFile.c_str());
    nrrdReader->Update();*/

    vtkSmartPointer<vtkMetaImageReader> metaReader = vtkSmartPointer<vtkMetaImageReader>::New();
    if (!metaReader->CanReadFile(srcFile.c_str()))
    {
        std::cerr << "Reader reports " << srcFile << " cannot be read.";
        return;
    }
    metaReader->SetFileName(srcFile.c_str());
    metaReader->Update();


    //1. read volume image with ITK whose InputPixelType is short
    using InputImageType = itk::Image<InputPixelType, 3>;
    using VTKToITKFilterType = itk::VTKImageToImageFilter<InputImageType>;
    auto filter = VTKToITKFilterType::New();
    filter->SetInput(metaReader->GetOutput());

    try
    {
        filter->Update();
    }
    catch (const itk::ExceptionObject & error)
    {
        std::cerr << "Error: " << error << std::endl;
        return;
    }
    InputImageType::ConstPointer itkInputImage = filter->GetOutput();

    using  OutputImageType = itk::Image<OutputPixelType, 3>;

    //2. using CastImageFilter to convert InputPixelType to OutputPixelType, short to float
    using CastFilterType = itk::CastImageFilter<InputImageType, OutputImageType>;
    CastFilterType::Pointer castFilter = CastFilterType::New();
    castFilter->SetInput(itkInputImage);
    try
    {
        castFilter->Update();
    }
    catch (itk::ExceptionObject e)
    {
        std::cout << e << std::endl;
    }

    itk::NrrdImageIO::Pointer io = itk::NrrdImageIO::New();
    io->SetFileType(itk::ImageIOBase::Binary);

    //set the nrrd meta data: space and space direction
    typedef itk::MetaDataDictionary DictionaryType;
    DictionaryType& dictionary = io->GetMetaDataDictionary();
    using MetaDataStringType = itk::MetaDataObject<std::string>;
    MetaDataStringType::Pointer spaceMeta = MetaDataStringType::New();
    spaceMeta->SetMetaDataObjectValue("scanner-xyz");
    MetaDataStringType::Pointer spaceMetaDirection = MetaDataStringType::New();
    spaceMetaDirection->SetMetaDataObjectValue("(0.25,0,0) (0,0.25,0) (0,0,1)");
    dictionary.Set("space", spaceMeta);
    dictionary.Set("space directions", spaceMetaDirection);
    


    typedef itk::ImageFileWriter<OutputImageType > WriterType;
    WriterType::Pointer nrrdWriter = WriterType::New();
    nrrdWriter->UseInputMetaDataDictionaryOn();
    nrrdWriter->SetInput(castFilter->GetOutput());
    nrrdWriter->SetImageIO(io);
    nrrdWriter->SetMetaDataDictionary(dictionary);
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


void addMetaTagsToNRRDFile(std::string srcFile, std::string outFile)
{
    ///////////////////////////////////////////////////////////////////////////////
    //if i don't put this line of code, it will report no registory for this kind type of file
    ///////////////////////////////////////////////////////////////////////////////
    itk::NrrdImageIOFactory::RegisterOneFactory();

    typedef itk::Image<float, 3> ImageType;
    typedef itk::ImageFileReader<ImageType> ReaderType;
    typedef itk::ImageFileWriter<ImageType> WriterType;
    typedef itk::MetaDataObject<std::string> MetaDataStringType;

    typedef itk::MetaDataObject<std::vector<double>> MetaDataDoubleArrayType;
    // add the meta tags
    //add space tag
    MetaDataStringType::Pointer meta = MetaDataStringType::New();
    meta->SetMetaDataObjectValue("scanner-xyz");
    //add the private tag
    MetaDataStringType::Pointer registrationType = MetaDataStringType::New();
    registrationType->SetMetaDataObjectValue("fixed");
    //add space direction
    MetaDataDoubleArrayType::Pointer spacingType = MetaDataDoubleArrayType::New();
    std::vector<double> directions = { 0.5, 0.5, 1,0.5, 0.5, 1,0.5, 0.5, 1 };
    spacingType->SetMetaDataObjectValue(directions);
    // directions
    using MatrixType = itk::Matrix<double, 3, 3>;
    typedef itk::MetaDataObject<MatrixType> MetaDataMatrixType;
    MatrixType M;
    M(0, 0) = 1.0;
    M(0, 1) = 2.0;
    M(0, 2) = 3.0;
    M(1, 0) = 4.0;
    M(1, 1) = 5.0;
    M(1, 2) = 6.0;
    M(2, 0) = 7.0;
    M(2, 1) = 8.0;
    M(2, 2) = 9.0;
    MetaDataMatrixType::Pointer directionsType = MetaDataMatrixType::New();
    directionsType->SetMetaDataObjectValue(M);


    ReaderType::Pointer reader = ReaderType::New();
    WriterType::Pointer writer = WriterType::New();
    reader->SetFileName(srcFile);
    reader->Update();


    using Dictionary = itk::MetaDataDictionary;
    using MetaDataStringType = itk::MetaDataObject<std::string>;
    Dictionary& dic = reader->GetOutput()->GetMetaDataDictionary();
    auto itr = dic.Begin();
    auto end = dic.End();
    //std::string entryId = "NRRD_space";
    std::string entryId = "ITK_original_spacing";
    auto tagItr = dic.Find(entryId);
    /*****************************************/
    /***********you can see all the meta data ****/
    /***********************************************/
    if (tagItr != end)
    {
        MetaDataStringType::ConstPointer entryvalue =
            dynamic_cast<const MetaDataStringType*> (tagItr->second.GetPointer());
        if (entryvalue)
        {
            std::string tagvalue = entryvalue->GetMetaDataObjectValue();
            std::cout << "patient: " << tagvalue << std::endl;
        }

        MetaDataDoubleArrayType::ConstPointer directionsValue =
            dynamic_cast<const MetaDataDoubleArrayType*> (tagItr->second.GetPointer());
        if (directionsValue)
        {
            std::vector<double> spacings = directionsValue->GetMetaDataObjectValue();
            std::cout << "spacings: " << spacings[0] << ", " << spacings[1] << ", " << spacings[2] << std::endl;
        }

        dic.Set("NRRD_space", meta);
        //dic.Set("registration type", registrationType);
        dic.Set("ITK_original_spacing", spacingType);

        dic.Set("ITK_original_direction", directionsType);
        MetaDataDoubleArrayType::ConstPointer directionsValue2 =
            dynamic_cast<const MetaDataDoubleArrayType*> (tagItr->second.GetPointer());
        if (directionsValue2)
        {
            std::vector<double> spacings = directionsValue2->GetMetaDataObjectValue();
            std::cout << "spacings: " << spacings[0] << ", " << spacings[1] << ", " << spacings[2] << std::endl;
        }
    }


    writer->SetFileName(outFile);
    writer->SetInput(reader->GetOutput());
    writer->Update();

    //this snippet of codes generate the nrrd file like this:
    /*
    NRRD0004
    # Complete NRRD file format specification at:
    # http://teem.sourceforge.net/nrrd/format.html
    type: float
    dimension: 3
    space: left-posterior-superior   //in itk, this space tag actually is NRRD_space
    sizes: 538 270 70
    space directions: (0.25,0,0) (0,0.25,0) (0,0,1)
    kinds: domain domain domain
    endian: little
    encoding: raw
    space origin: (67.233900000000006,67.490700000000004,0)
    ITK_InputFilterName:=NrrdImageIO
    ITK_original_direction:=
    ITK_original_spacing:=
    space:=scanner-xyz
    */
    // the tag is built like that, i don't think it generate the 

}