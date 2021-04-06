/*=========================================================================
 *
 *  Copyright SINAPSE: Scalable Informatics for Neuroscience, Processing and Software Engineering
 *            The University of Iowa
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
#include <sstream>
#include "itkMedianImageFilter.h"
#include "itkExtractImageFilter.h"
#include "BRAINSCommonLib.h"
#include "BRAINSThreadControl.h"
#include "BRAINSFitHelper.h"

#include "BRAINSToolsVersion.h"


// This program was modified from
// Insight/Examples/Registration/ImageRegistration8.cxx
// and is an improved replacement for the old (and defective)

using BRAINSFitPixelType = float;
using FixedVolumeType = itk::Image<BRAINSFitPixelType, Dimension>;
using MovingVolumeType = itk::Image<BRAINSFitPixelType, Dimension>;

using InputImageType = itk::Image<BRAINSFitPixelType, MaxInputDimension>;
using FixedVolumeReaderType = itk::ImageFileReader<InputImageType>;
using MovingVolumeReaderType = itk::ImageFileReader<InputImageType>;
using AffineTransformPointer = itk::AffineTransform<double, 3>::Pointer;

template <typename ImageType>
typename ImageType::Pointer
ExtractImage(typename InputImageType::Pointer & inputImage, unsigned int InputImageTimeIndex)
{
  using ExtractImageFilterType = typename itk::ExtractImageFilter<InputImageType, ImageType>;
  typename ExtractImageFilterType::Pointer extractImageFilter = ExtractImageFilterType::New();
  extractImageFilter->SetDirectionCollapseToSubmatrix();

  // fixedVolumeReader->GetOutput();
  InputImageType::RegionType inputRegion = inputImage->GetLargestPossibleRegion();
  InputImageType::SizeType   inputSize = inputRegion.GetSize();
  inputSize[3] = 0;
  inputRegion.SetSize(inputSize);

  InputImageType::IndexType inputIndex = inputRegion.GetIndex();
  inputIndex[0] = 0;
  inputIndex[1] = 0;
  inputIndex[2] = 0;
  inputIndex[3] = InputImageTimeIndex;
  inputRegion.SetIndex(inputIndex);
  extractImageFilter->SetExtractionRegion(inputRegion);
  extractImageFilter->SetInput(inputImage);

  try
  {
    extractImageFilter->Update();
  }
  catch (...)
  {
    std::cout << "Error while extracting a time indexed fixed image." << std::endl;
    throw;
  }

  typename ImageType::Pointer extractImage = extractImageFilter->GetOutput();
  //  std::cerr << "Extract fixed image origin" << extractImage->GetOrigin() << std::endl;

  return extractImage;
}

template <typename ImageType>
typename ImageType::Pointer
DoMedian(typename ImageType::Pointer & input, typename ImageType::SizeType indexRadius)
{
  using MedianFilterType = typename itk::MedianImageFilter<ImageType, ImageType>;
  typename MedianFilterType::Pointer medianFilter = MedianFilterType::New();

  medianFilter->SetRadius(indexRadius);
  medianFilter->SetInput(input);
  medianFilter->Update();
  typename ImageType::Pointer result = medianFilter->GetOutput();
  return result;
}

#ifdef USE_DebugImageViewer
/*************************
 * Have a global variable to
 * add debugging information.
 */
DebugImageViewerClient DebugImageDisplaySender;
#endif

int
main(int argc, char * argv[])
{
	int numberOfThreads = -1;
	using GenericTransformType = itk::Transform<double, 3, 3>;
	const BRAINSUtils::StackPushITKDefaultNumberOfThreads TempDefaultNumberOfThreadsHolder(numberOfThreads);

	std::string fixedVolume = argv[1];
	std::string movingVolume = argv[2];
	// Extracting a timeIndex cube from the fixed image goes here....
	// Also MedianFilter
	FixedVolumeType::Pointer  extractFixedVolume;
	MovingVolumeType::Pointer extractMovingVolume;
	unsigned int fixedVolumeTimeIndex = 0;
	unsigned int movingVolumeTimeIndex = 0;
	InputImageType::Pointer   OriginalFixedVolume(itkUtil::ReadImage<InputImageType>(fixedVolume));

	std::cout << "Original Fixed image origin" << OriginalFixedVolume->GetOrigin() << std::endl;
	// fixedVolumeTimeIndex lets lets us treat 3D as 4D.
	/***********************
	* Acquire Fixed Image Index
	**********************/
	extractFixedVolume = ExtractImage<FixedVolumeType>(OriginalFixedVolume, fixedVolumeTimeIndex);
	// Extracting a timeIndex cube from the moving image goes here....

	InputImageType::Pointer OriginalMovingVolume(itkUtil::ReadImage<InputImageType>(movingVolume));
	// This default lets us treat 3D as 4D.
	// const unsigned int movingVolumeTimeIndex;

	/***********************
	* Acquire Moving Image Index
	**********************/
	extractMovingVolume = ExtractImage<MovingVolumeType>(OriginalMovingVolume, movingVolumeTimeIndex);



	using HelperType = itk::BRAINSFitHelper;
	HelperType::Pointer myHelper = HelperType::New();

	std::vector<std::string> localTransformType;
	localTransformType.emplace_back("Rigid");
	myHelper->SetTransformType(localTransformType);

	myHelper->SetFixedVolume(extractFixedVolume);
	myHelper->SetMovingVolume(extractMovingVolume);
	bool histogramMatch = false;
	myHelper->SetHistogramMatch(histogramMatch);
	int removeIntensityOutliers = 0;
	myHelper->SetRemoveIntensityOutliers(removeIntensityOutliers);
	int numberOfMatchPoints = 10;
	myHelper->SetNumberOfMatchPoints(numberOfMatchPoints);
	//myHelper->SetFixedBinaryVolume(fixedMask);
	//myHelper->SetMovingBinaryVolume(movingMask);
	//myHelper->SetOutputFixedVolumeROI(outputFixedVolumeROI);
	//myHelper->SetOutputMovingVolumeROI(outputMovingVolumeROI);
	double samplingPercentage = 0.002;
	myHelper->SetSamplingPercentage(samplingPercentage);
	double numberOfHistogramBins = 50;
	myHelper->SetNumberOfHistogramBins(numberOfHistogramBins);
	std::vector<int> numberOfIterations(1);
	numberOfIterations[0] = 1500;
	myHelper->SetNumberOfIterations(numberOfIterations);
	//std::vector<double> maximumStepLength(4);
	//maximumStepLength[0] = 0.05;
	double maximumStepLength = 0.05;
	myHelper->SetMaximumStepLength(maximumStepLength);
	std::vector<double> minimumStepLength(1);
	minimumStepLength[0] = 0.001;
	//double minimumStepLength = 0.001;
	myHelper->SetMinimumStepLength(minimumStepLength);
	double relaxationFactor = 0.5;
	myHelper->SetRelaxationFactor(relaxationFactor);
	double translationScale = 1000.0;
	myHelper->SetTranslationScale(translationScale);
	double reproportionScale = 1.0;
	myHelper->SetReproportionScale(reproportionScale);
	double skewScale = 1.0;
	myHelper->SetSkewScale(skewScale);
	//myHelper->SetBackgroundFillValue(backgroundFillValue);
	std::string  localInitializeTransformMode = "Off";
	myHelper->SetInitializeTransformMode(localInitializeTransformMode);
	double maskInferiorCutOffFromCenter = 1000;
	myHelper->SetMaskInferiorCutOffFromCenter(maskInferiorCutOffFromCenter);
	using CompositeTransformType = itk::CompositeTransform<double, 3>;
	CompositeTransformType::Pointer currentGenericTransform = nullptr;
	myHelper->SetCurrentGenericTransform(currentGenericTransform);
	//myHelper->SetSplineGridSize(BSplineGridSize);
	//myHelper->SetCostFunctionConvergenceFactor(costFunctionConvergenceFactor);
	//myHelper->SetProjectedGradientTolerance(projectedGradientTolerance);
	//myHelper->SetMaxBSplineDisplacement(maxBSplineDisplacement);
	bool UseDebugImageViewer = true;
	myHelper->SetDisplayDeformedImage(UseDebugImageViewer);
	//myHelper->SetPromptUserAfterDisplay(PromptAfterImageSend);
	int debugLevel = 0;
	myHelper->SetDebugLevel(debugLevel);
	std::string costMetric = "MMI";	
	//std::string costMetric = "MSE";
	myHelper->SetCostMetricName(costMetric);
	int useROIBSpline = false;
	myHelper->SetUseROIBSpline(useROIBSpline);
	std::string metricSamplingStrategy = "Random";
	myHelper->SetSamplingStrategy(metricSamplingStrategy);
	bool initializeRegistrationByCurrentGenericTransform = false;
	myHelper->SetInitializeRegistrationByCurrentGenericTransform(initializeRegistrationByCurrentGenericTransform); 
	int maximumNumberOfEvaluations = 900;
	myHelper->SetMaximumNumberOfEvaluations(maximumNumberOfEvaluations);
	int maximumNumberOfCorrections = 25;
	myHelper->SetMaximumNumberOfCorrections(maximumNumberOfCorrections);
	std::vector<int> BSplineGridSize{ 14,10,12 };
	myHelper->SetSplineGridSize(BSplineGridSize);
	bool writeOutputTransformInFloat = false;
	myHelper->SetWriteOutputTransformInFloat(writeOutputTransformInFloat);

	try
	{
		myHelper->Update();
	}
	catch (itk::ExceptionObject & err)
	{
		std::cerr << "Exception during registration: " << std::endl;
		std::cerr << err << std::endl;
		return EXIT_FAILURE;
	}
	currentGenericTransform = myHelper->GetCurrentGenericTransform();
	//const itk::Transform<double, 3, 3>::ConstPointer currInitTransformFormGenericComposite =
	//	currentGenericTransform->GetFrontTransform();
	//std::string type = currInitTransformFormGenericComposite->GetTransformTypeAsString();
	//std::cout << "type: " << type << std::endl;
	//currInitTransformFormGenericComposite->GetParameters();

	using TransformType = itk::VersorRigid3DTransform<double>;
	TransformType::Pointer finalTransform = TransformType::New();

	finalTransform->SetFixedParameters(
		currentGenericTransform->GetFixedParameters());
	finalTransform->SetParameters(currentGenericTransform->GetParameters());

	// Software Guide : BeginCodeSnippet
	TransformType::MatrixType matrix = finalTransform->GetMatrix();
	TransformType::OffsetType offset = finalTransform->GetOffset();
	TransformType::TranslationType translation = finalTransform->GetTranslation();
	std::cout << "Matrix = " << std::endl << matrix << std::endl;
	std::cout << "Offset = " << std::endl << offset << std::endl;
	std::cout << "translation = " << std::endl << translation << std::endl;

	std::cout << "hello world!" << std::endl;
}
