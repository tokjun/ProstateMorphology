#include "itkImageFileWriter.h"

#include "itkSmoothingRecursiveGaussianImageFilter.h"
#include "itkLabelStatisticsImageFilter.h"
#include "itkSignedMaurerDistanceMapImageFilter.h"
#include "itkConnectedComponentImageFilter.h"
#include "itkRelabelComponentImageFilter.h"
#include "itkCastImageFilter.h"

#include "itkPluginUtilities.h"

#include "ProstateMorphologyCLP.h"

// Use an anonymous namespace to keep class types and function names
// from colliding when module is used as shared object module.  Every
// thing should be in an anonymous namespace except for the module
// entry point, e.g. main()
//
namespace
{

template <class T>
int DoIt( int argc, char * argv[], T )
{
  PARSE_ARGS;

  typedef    T InputPixelType;
  typedef    T OutputPixelType;
  typedef    unsigned int LabelInputPixelType;
  typedef    float DistanceMapPixelType;

  typedef itk::Image<InputPixelType,  3> InputImageType;
  typedef itk::Image<OutputPixelType, 3> OutputImageType;
  typedef itk::Image<DistanceMapPixelType, 3> DistanceMapImageType;

  typedef itk::Image<LabelInputPixelType, 3> LabelImageType;
  
  typedef itk::ImageFileReader<InputImageType>  ReaderType;
  typedef itk::ImageFileWriter<OutputImageType> WriterType;
  typedef itk::ImageFileReader<LabelImageType>  LabelReaderType;
  
  typedef itk::LabelStatisticsImageFilter<
    InputImageType, LabelImageType > LabelStatisticsImageFilterType;
  
  typedef itk::LabelStatisticsImageFilter<
    DistanceMapImageType, LabelImageType > DistanceMapLabelStatisticsImageFilterType;
  
  typedef itk::SmoothingRecursiveGaussianImageFilter<
    InputImageType, OutputImageType>  FilterType;

  typedef   itk::ConnectedComponentImageFilter<
    LabelImageType, LabelImageType >  CCFilterType;
  
  typedef   itk::RelabelComponentImageFilter<
    LabelImageType, LabelImageType > RelabelType;
  
  typedef   itk::SignedMaurerDistanceMapImageFilter<
    LabelImageType, DistanceMapImageType >  DistanceMapFilterType;

  typedef itk::CastImageFilter< LabelImageType, OutputImageType > CastToOutputFilterType;  

  // Configure input image reader
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( inputVolume.c_str() );

  // Configure input label readers
  typename LabelReaderType::Pointer tgReader = LabelReaderType::New();
  tgReader->SetFileName( TG.c_str() );

  typename LabelReaderType::Pointer nvbReader = LabelReaderType::New();
  nvbReader->SetFileName( NVB.c_str() );

  typename LabelReaderType::Pointer eusReader = LabelReaderType::New();
  eusReader->SetFileName( EUS.c_str() );

  typename LabelReaderType::Pointer svReader = LabelReaderType::New();
  svReader->SetFileName( SV.c_str() );

  typename LabelReaderType::Pointer tumorReader = LabelReaderType::New();
  tumorReader->SetFileName( Tumor.c_str() );


  typename LabelStatisticsImageFilterType::Pointer lsiFilter = LabelStatisticsImageFilterType::New();
  lsiFilter->SetInput( reader->GetOutput() );

  // Run label statistics filter on each label
  typedef std::vector < LabelReaderType::Pointer > ReaderListType;

  // Re-label NVB
  typename CCFilterType::Pointer CCFilter = CCFilterType::New();
  typename RelabelType::Pointer relabelFilter = RelabelType::New();
  CCFilter->SetInput( nvbReader->GetOutput() );
  CCFilter->FullyConnectedOff();
  relabelFilter->SetInput ( CCFilter->GetOutput() );  
  //relabelFilter->SetMinimumObjectSize(2); // Remove noise
  relabelFilter->Update();

  // Generate a tumor distance map
  typename DistanceMapFilterType::Pointer distanceFilter = DistanceMapFilterType::New();
  distanceFilter->SetInput(tumorReader->GetOutput());
  distanceFilter->SquaredDistanceOff();
  distanceFilter->SetUseImageSpacing(true);
  try
    {
    distanceFilter->Update();
    }
  catch (itk::ExceptionObject &err)
    {
    std::cerr << err << std::endl;
    return EXIT_FAILURE ;
    }

  typename DistanceMapLabelStatisticsImageFilterType::Pointer dlsiFilter
    = DistanceMapLabelStatisticsImageFilterType::New();
  dlsiFilter->SetInput( distanceFilter->GetOutput() );
  dlsiFilter->SetLabelInput( relabelFilter->GetOutput() );
  dlsiFilter->Update();

  //typedef typename DistanceMapLabelStatisticsImageFilterType::ValidLabelValuesContainerType ValidLabelValuesType;
  //typedef typename DistanceMapLabelStatisticsImageFilterType::LabelPixelType                LabelPixelType;
  //typename ValidLabelValuesType::const_iterator vIt;
  //std::cout << "Minimum distance between tumor and NVB:" << std::endl;
  //for(vIt=dlsiFilter->GetValidLabelValues().begin();
  //    vIt != dlsiFilter->GetValidLabelValues().end();
  //    ++vIt)
  //  {
  //  if ( *vIt!= 0 && dlsiFilter->HasLabel(*vIt) )
  //    {
  //    LabelPixelType labelValue = *vIt;
  //    std::cout << "label: " << labelValue << std::endl;
  //    std::cout << "min: " << dlsiFilter->GetMinimum( labelValue ) << std::endl;
  //    //std::cout << "max: " << lsiFilter->GetMaximum( labelValue ) << std::endl;
  //    //std::cout << "median: " << lsiFilter->GetMedian( labelValue ) << std::endl;
  //    //std::cout << "mean: " << lsiFilter->GetMean( labelValue ) << std::endl;
  //    //std::cout << "sigma: " << lsiFilter->GetSigma( labelValue ) << std::endl;
  //    //std::cout << "variance: " << lsiFilter->GetVariance( labelValue ) << std::endl;
  //    //std::cout << "sum: " << lsiFilter->GetSum( labelValue ) << std::endl;
  //    std::cout << "count: " << lsiFilter->GetCount( labelValue ) << std::endl;
  //    //std::cout << "region: " << lsiFilter->GetRegion( labelValue ) << std::endl;
  //    std::cout << std::endl;
  //    }
  //  }

  // Measure Volumes
  float volTG;
  float volEUS;
  float volTumor;
  float volSV;
  float distTG;
  float distEUS;
  float distSV;

  const int labelTG = 1;
  const int labelSV = 1;
  const int labelEUS = 1;
  const int labelTumor = 1;
  const int labelNVB1 = 1;
  const int labelNVB2 = 2;
  const int labelBackground = 0;


  typename InputImageType::SpacingType spacing;

  // Check if there are two NVBs
  float distNVB[2];
  if (dlsiFilter->HasLabel(labelNVB1) && dlsiFilter->HasLabel(labelNVB2))
    {
    distNVB[0] = dlsiFilter->GetMinimum( labelNVB1 );
    distNVB[1] = dlsiFilter->GetMinimum( labelNVB2 );
    }
  else
    {
    std::cerr << "Failed to detect NVBs on both side." << std::endl;
    }
  
  // TG
  dlsiFilter->SetLabelInput(tgReader->GetOutput());
  dlsiFilter->Update();
  if (dlsiFilter->HasLabel(labelTG))
    {
    spacing = tgReader->GetOutput()->GetSpacing();
    volTG = dlsiFilter->GetCount(labelTG) * spacing[0] * spacing[1] * spacing[2] / 1000.0;
    distTG = dlsiFilter->GetMinimum(labelBackground); // distance in background to check the distance from the surface
    }
  else
    {
    std::cerr << "Failed to detect TG." << std::endl;
    }
  
  // EUS
  dlsiFilter->SetLabelInput(eusReader->GetOutput());
  dlsiFilter->Update();
  if (dlsiFilter->HasLabel(labelEUS))
    {
    spacing = eusReader->GetOutput()->GetSpacing();
    volEUS = dlsiFilter->GetCount(labelEUS) * spacing[0] * spacing[1] * spacing[2] / 1000.0;
    distEUS = dlsiFilter->GetMinimum( labelEUS );
    }
  else
    {
    std::cerr << "Failed to detect EUS." << std::endl;
    }

  // SV
  dlsiFilter->SetLabelInput(svReader->GetOutput());
  dlsiFilter->Update();
  if (dlsiFilter->HasLabel(labelSV))
    {
    spacing = svReader->GetOutput()->GetSpacing();
    volSV = dlsiFilter->GetCount(labelSV) * spacing[0] * spacing[1] * spacing[2] / 1000.0;
    distSV = dlsiFilter->GetMinimum( labelSV );
    }
  else
    {
    std::cerr << "Failed to detect SV." << std::endl;
    }

  // Tumor
  dlsiFilter->SetLabelInput(tumorReader->GetOutput());
  dlsiFilter->Update();
  if (dlsiFilter->HasLabel(labelTumor))
    {
    spacing = tumorReader->GetOutput()->GetSpacing();
    volTumor = dlsiFilter->GetCount(labelTumor) * spacing[0] * spacing[1] * spacing[2]  / 1000.0;
    }
  else
    {
    std::cerr << "Failed to detect Tumor." << std::endl;
    }

  std::cout << "(volTG/volEUS/volSV/volTumor/distNVB1/distNVB2/distTG/distEUS/distSV),"
            << volTG << ", " << volEUS << ", " << volSV << ", " << volTumor << ", "
            << distNVB[0] << ", " << distNVB[1] << ","
            << distTG << ", " << distEUS << "," << distSV << std::endl;

  typename CastToOutputFilterType::Pointer castFilter = CastToOutputFilterType::New();
  castFilter->SetInput( relabelFilter->GetOutput() );

  typename WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( outputVolume.c_str() );
  writer->SetInput( castFilter->GetOutput() );
  writer->SetUseCompression(1);
  writer->Update();

  return EXIT_SUCCESS;
}

} // end of anonymous namespace

int main( int argc, char * argv[] )
{
  PARSE_ARGS;

  itk::ImageIOBase::IOPixelType     pixelType;
  itk::ImageIOBase::IOComponentType componentType;

  try
    {
    itk::GetImageType(inputVolume, pixelType, componentType);

    // This filter handles all types on input, but only produces
    // signed types
    switch( componentType )
      {
      case itk::ImageIOBase::UCHAR:
        return DoIt( argc, argv, static_cast<unsigned char>(0) );
        break;
      case itk::ImageIOBase::CHAR:
        return DoIt( argc, argv, static_cast<char>(0) );
        break;
      case itk::ImageIOBase::USHORT:
        return DoIt( argc, argv, static_cast<unsigned short>(0) );
        break;
      case itk::ImageIOBase::SHORT:
        return DoIt( argc, argv, static_cast<short>(0) );
        break;
      case itk::ImageIOBase::UINT:
        return DoIt( argc, argv, static_cast<unsigned int>(0) );
        break;
      case itk::ImageIOBase::INT:
        return DoIt( argc, argv, static_cast<int>(0) );
        break;
      case itk::ImageIOBase::ULONG:
        return DoIt( argc, argv, static_cast<unsigned long>(0) );
        break;
      case itk::ImageIOBase::LONG:
        return DoIt( argc, argv, static_cast<long>(0) );
        break;
      case itk::ImageIOBase::FLOAT:
        return DoIt( argc, argv, static_cast<float>(0) );
        break;
      case itk::ImageIOBase::DOUBLE:
        return DoIt( argc, argv, static_cast<double>(0) );
        break;
      case itk::ImageIOBase::UNKNOWNCOMPONENTTYPE:
      default:
        std::cout << "unknown component type" << std::endl;
        break;
      }
    }

  catch( itk::ExceptionObject & excep )
    {
    std::cerr << argv[0] << ": exception caught !" << std::endl;
    std::cerr << excep << std::endl;
    return EXIT_FAILURE;
    }
  return EXIT_SUCCESS;
}
