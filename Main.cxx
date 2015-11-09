#include "itkCastImageFilter.h"

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include <iostream>
#include <sstream>
#include <string>


#include <metaCommand.h>

#include "itkImage.h"

#include <errno.h>
#include <iostream>
#include <limits.h>
#include "itkExponentialDeformationFieldImageFilter2.h"
#include "itksys/SystemTools.hxx"
#include "itkRegularExpressionSeriesFileNames.h"
#include "itkNumericSeriesFileNames.h"
#include "itkExtractImageFilter.h"

#include "rpiLCClogDemons.hxx"

#include "itksys/SystemTools.hxx"

struct arguments
{

    std::string inputImage;                   /* -I option */
    std::string maskImage;                     /* -M option */
    int nbrIterations;                       /* -t option */
    int updateRule;                         /* -r option */
    int regularization;                     /* -R option */
    double sigmaI;                          /* -S option */
    double sigmaVel;                        /* -d option */
    std::string outputFolder;                   /* -o option */
    int refMesh;                                /* -refMesh option*/
    bool alphabeticalSort;                      /* -alphabeticalSort option*/
    std::string regularExpression;              /*-reg option*/


    friend std::ostream& operator<< (std::ostream& o, const arguments& args)
    {

        return o
        <<"  Arguments structure:"<<std::endl
        <<" Mask image path:" << std::endl
        <<" Number Iterations: "<<args.nbrIterations<<std::endl
        <<" Update Rule:" << args.updateRule<< std::endl
        <<" Regularization:" << args.regularization<< std::endl
        <<" Sigma I:" << args.sigmaI<< std::endl
        <<" Sigma Vel:" << args.sigmaVel<< std::endl
        <<" Input Folder:" << args.inputImage << std::endl
        <<" Reg:" << args.regularExpression << std::endl
        <<" Ref Mesh:" << args.refMesh << std::endl
        <<" Alphabetical Sort:" << args.alphabeticalSort << std::endl
        <<" Output Folder:" << args.outputFolder << std::endl;

    }
};

void help_callback()
{
    std::cout<<std::endl;
    std::cout<<"Copyright (c) 2015 INRIA"<<std::endl;
    std::cout<<"Code: Marc-Michel Rohe"<<std::endl;
    std::cout<<"Report bugs to <marc-michel.rohe \\at inria.fr>"<<std::endl;

    exit( EXIT_FAILURE );
}

int atoi_check( const char * str )
{
  char *endptr;
  long val= strtol(str, &endptr, 0);

  /* Check for various possible errors */
  if ( (errno == ERANGE && (val == LONG_MAX || val == LONG_MIN))
     || val>=INT_MAX || val<=INT_MIN )
    {
    std::cout<<std::endl;
    std::cout<<"Cannot parse integer. Out of bound."<<std::endl;
    exit( EXIT_FAILURE );
    }

  if (endptr == str || *endptr!='\0')
    {
    std::cout<<std::endl;
    std::cout<<"Cannot parse integer. Contains non-digits or is empty."<<std::endl;
    exit( EXIT_FAILURE );
    }

  return val;
}

std::vector<unsigned int> parseUIntVector( const std::string & str)
{
  std::vector<unsigned int> vect;

  std::string::size_type crosspos = str.find('x',0);

  if (crosspos == std::string::npos)
    {
    // only one uint
    vect.push_back( static_cast<unsigned int>( atoi_check(str.c_str()) ));
    return vect;
    }

  // first uint
  vect.push_back( static_cast<unsigned int>(
                    atoi_check( (str.substr(0,crosspos)).c_str()  ) ));

  while(true)
    {
    std::string::size_type crossposfrom = crosspos;
    crosspos =  str.find('x',crossposfrom+1);

    if (crosspos == std::string::npos)
      {
      vect.push_back( static_cast<unsigned int>(
                         atoi_check( (str.substr(crossposfrom+1,str.length()-crossposfrom-1)).c_str()  ) ));
      return vect;
      }

    vect.push_back( static_cast<unsigned int>(
                       atoi_check( (str.substr(crossposfrom+1,crosspos-crossposfrom-1)).c_str()  ) ));
    }
}

void parseOpts (int argc, char **argv, struct arguments & args)
{
    // Command line parser
    MetaCommand command;
    command.SetParseFailureOnUnrecognizedOption( true );
    command.SetHelpCallBack(help_callback);



    // Define parsing options



    command.SetOption("Iteration","t",false,"Iteration");
    command.SetOptionLongTag("Iteration","nbr-Iteration");
    command.AddOptionField("Iteration","intval",MetaCommand::INT,false,"3");

    command.SetOption("OutputFolder","o",false,"Ouput folder");
    command.SetOptionLongTag("OutputFolder","output-folder");
    command.AddOptionField("OutputFolder","filename",MetaCommand::STRING,false,"");

    command.SetOption("InputImage","I",true,"Input folder");
    command.SetOptionLongTag("InputImage","input-folder");
    command.AddOptionField("InputImage","filename",MetaCommand::STRING,true,"");

    command.SetOption("UpdateRule","r",false,"Update rule:  0 : SSD similarity - exp(v)<-exp(v)oexp(u) (log-domain),1 : SSD similarity - exp(v)<-symmetrized(exp(v)oexp(u)) (symmetriclog-domain), 2 : LCC similarity (default 1).");
    command.SetOptionLongTag("UpdateRule","update-rule");
    command.AddOptionField("UpdateRule","intval",MetaCommand::INT,false,"1");

    command.SetOption("regularization","R",false,"Type of Regularization: 0 = Gaussian convolution (classical Demons), 1 = Harmonic + Bending Energy (default 1)");
    command.SetOptionLongTag("regularization","type-reg");
    command.AddOptionField("regularization","intval",MetaCommand::INT,false,"0");

    command.SetOption("maskImage","M",false,"Path of the mask image");
    command.SetOptionLongTag("maskImage","mask-image");
    command.AddOptionField("maskImage","filename",MetaCommand::STRING,false);

    command.SetOption("SigmaI","S",false,"Trade-off between similarity and regularization (sigma_i ^2): 0 (sharper but unregular deformations) < sigma_i <= 1 (smootherdeformations but weaker correspondencies). Default: 0.15");
    command.SetOptionLongTag("SigmaI","sigma-I");
    command.AddOptionField("SigmaI","floatval",MetaCommand::FLOAT,false,"0.5");

    command.SetOption("SigmaVel","d",false,"Standard deviation of the Gaussian smoothing of the stationary velocity field (world units). Setting it below 0.1 means no smoothing will be performed (default 1.5).");
    command.SetOptionLongTag("SigmaVel","sigma-vel");
    command.AddOptionField("SigmaVel","floatval",MetaCommand::FLOAT,false,"4");

    command.SetOption("Reg","",false,"Standard deviation of the Gaussian smoothing of the stationary velocity field (world units). Setting it below 0.1 means no smoothing will be performed (default 1.5).");
    command.SetOptionLongTag("Reg","reg");
    command.AddOptionField("Reg","filename",MetaCommand::FLOAT,false,"([0-9]+)[.]");

    command.SetOption("RefMesh","",false,"Ref Mesh of the series generator");
    command.SetOptionLongTag("RefMesh","refMesh");
    command.AddOptionField("RefMesh","intval",MetaCommand::INT,false,"1");

    command.SetOption("AlphabeticalSort","",false,"Alphabetical Sort of the Series Generatoc ");
    command.SetOptionLongTag("AlphabeticalSort","alphabeticalSort");
    command.AddOptionField("AlphabeticalSort","boolval",MetaCommand::BOOL,false,"0");

    // Actually parse the command line
    if (!command.Parse(argc,argv))
    {
        exit( EXIT_FAILURE );
    }


    // Store the parsed information into a struct

    args.outputFolder = command.GetValueAsString("OutputFolder","filename");
    args.inputImage = command.GetValueAsString("InputImage","filename");
    args.nbrIterations = command.GetValueAsInt("Iteration","intval");
    args.updateRule= command.GetValueAsInt("UpdateRule","intval");
    args.regularization= command.GetValueAsInt("regularization","intval");
    args.sigmaI= command.GetValueAsFloat("SigmaI","floatval");
    args.sigmaVel= command.GetValueAsFloat("SigmaVel","floatval");
    args.refMesh=command.GetValueAsInt("RegMesh","intval");
            args.regularExpression=command.GetValueAsFloat("Reg","filename");
            args.alphabeticalSort=command.GetValueAsBool("AlphabeticalSort","boolval");
            args.maskImage=command.GetValueAsString("maskImage","filename");

}




int main( int argc, char *argv[] )
{
    typedef float TScalar;
    const int Dimension = 3;
    typedef itk::Vector<double,Dimension>                    VectorType;
    typedef itk::Vector<double,Dimension-1>                    SliceVectorType;

    typedef itk::Image<VectorType,Dimension>                FieldType;
    typedef itk::Image<SliceVectorType,Dimension-1>                SliceFieldType;
    typedef itk::ImageFileReader<FieldType>                FieldTypeReader;
    typedef itk::ImageFileWriter<FieldType>                FieldTypeWriter;


    typedef itk::Image<double,Dimension>                   ImageType;
    typedef itk::Image<double, Dimension-1>                 SliceImageType;

    typedef ImageType::Pointer                              ImagePointer;
    typedef itk::ImageFileWriter<ImageType>                   ImageFileWriter;
        typedef ImageFileWriter::Pointer                              ImageFileWriterPointer;

    typedef itk::ImageFileReader<ImageType>       ImageReaderType;
    typedef ImageReaderType::Pointer                        ImageReaderPointer;

    typedef itk::ExtractImageFilter<ImageType,SliceImageType> ExtractSliceFilterType;

    typedef itk::StationaryVelocityFieldTransform<double, 3>  FieldTransformType;
    typedef itk::StationaryVelocityFieldTransform<double, 2>  SliceFieldTransformType;


    typedef itk::ExponentialDeformationFieldImageFilter<SliceFieldType,SliceFieldType> ExponentialFilterType;
    typedef ExponentialFilterType::Pointer                                   ExponentialFilterPointer;

    typedef itk::WarpImageFilter< SliceImageType, SliceImageType, SliceFieldType >  WarperType;

    struct arguments args;
    parseOpts (argc, argv, args);

    /// Vector type
    typedef vnl_vector<TScalar> VNLVectorType;
    /// Matrix type.
    typedef vnl_matrix<TScalar> VNLMatrixType;
    /// List of matrices type.
    typedef std::vector<VNLMatrixType> VNLMatrixList;

    std::cout << "Starting Resample Barycentric with the following arguments:" << std::endl;
    std::cout<<args<<std::endl<<std::endl;


    /*

    typedef itk::RegularExpressionSeriesFileNames    NameGeneratorType;
    NameGeneratorType::Pointer nameGenerator = NameGeneratorType::New();

    std::string input_directory = args.inputFolder;
    std::string regularExpression = args.regularExpression;
    unsigned int subMatch = args.refMesh;
    if (args.alphabeticalSort)
        nameGenerator->NumericSortOn();
    else
        nameGenerator->NumericSortOff();

    //  Prep. the generator
    nameGenerator->SetRegularExpression( regularExpression );
    nameGenerator->SetSubMatch( subMatch );
    nameGenerator->SetDirectory( input_directory );

    // We are ready to read all the images. Print the filenames too.


    unsigned int numberOfFilenames =  nameGenerator->GetFileNames().size();

    for(unsigned int fni = 0; fni<numberOfFilenames; fni++)
    {
        std::cout << "filename # " << fni << " = ";
        std::cout << nameGenerator->GetFileNames()[fni] << std::endl;
    }

    // We can trigger the reading process by calling the \code{Update()} method on
    // the reader. It is wise to put this invocation inside a
    // \code{try/catch} block since the process may eventually throw exceptions.




    //  Now we can prepare the process for writing the dataset. First, we take the
    //  name of the output directory from the command line arguments.

    std::string output_directory = args.outputFolder;

    //  Second, we make sure the output directory exist, using the cross platform
    //  tools: itksys::SystemTools. In this case we select to create the directory
    //  if it does not exist yet.
    //
    //  \index{itksys!SystemTools}
    //  \index{itksys!MakeDirectory}
    //  \index{SystemTools}
    //  \index{SystemTools!MakeDirectory}
    //  \index{MakeDirectory!SystemTools}
    //  \index{MakeDirectory!itksys}

    itksys::SystemTools::MakeDirectory( output_directory.c_str() );


    typedef itk::NumericSeriesFileNames				OutputNamesGeneratorType;

    std::string format = args.outputFolder + "/image%03d.mha";

    OutputNamesGeneratorType::Pointer outputNamesGenerator = OutputNamesGeneratorType::New();
    outputNamesGenerator->SetSeriesFormat( format.c_str() );
    outputNamesGenerator->SetStartIndex( 0 );
    outputNamesGenerator->SetEndIndex( numberOfFilenames-1 );
    outputNamesGenerator->SetIncrementIndex( 1 );


    */

    ImagePointer maskImage =0;


    if (!args.maskImage.empty())
    {
        ImageReaderPointer maskReader=ImageReaderType::New();
        maskReader->SetFileName(args.maskImage);
        maskReader->Update();

        maskImage = maskReader->GetOutput();
        maskImage->DisconnectPipeline();
    }





    //for(unsigned int fni = 0; fni<numberOfFilenames; fni++)
    {
        std::cout << "<<<<<<====================>>>> " << std::endl;
        std::cout << "<<<<<<====================>>>> " << std::endl;
        std::cout << "<<<<<<====================>>>> " << std::endl;

        //std::cout << "Processing Image number: " << fni << std::endl;
        ImageReaderType::Pointer reader = ImageReaderType::New();
        //reader->SetFileName( nameGenerator->GetFileNames()[fni] );
        reader->SetFileName( args.inputImage);


        try { reader->Update(); }
        catch (itk::ExceptionObject &excp) {
            std::cerr << "Exception thrown while reading the images" << std::endl;
            std::cerr << excp << std::endl;
            return EXIT_FAILURE;
        }

        ImageType::Pointer image = ImageType::New();
        image=reader->GetOutput();
        image->DisconnectPipeline();


        ImageType::Pointer resampledImage = ImageType::New();
        resampledImage->CopyInformation(image);

        ImageType::SpacingType baseSpacing;
        baseSpacing=image->GetSpacing();

        ImageType::SizeType size;
        size=image->GetLargestPossibleRegion().GetSize();
        int sizez= (int) (baseSpacing[2]/(double)baseSpacing[0])*(size[2]-1)+1;

        sizez=sizez-1;
        size[2]=sizez;

        ImageType::SpacingType spacing;
        spacing=baseSpacing;
        spacing[2]=baseSpacing[0];
        resampledImage->SetSpacing(spacing);

        ImageType::IndexType index;
        index.Fill(0);

        ImageType::RegionType region;
        region.SetSize(size);
        region.SetIndex(index);

        resampledImage->SetRegions(region);
        resampledImage->Allocate();

        const SliceFieldType * svf2= SliceFieldType::New();
        const SliceFieldType * svf1= SliceFieldType::New();

        for (int j=0; j<sizez ; ++j)
        {


            std::cout << "Lambda registration for slice " << j << std::endl;
            int sliceDown= (int) j*spacing[2]/baseSpacing[2];
            int sliceUp=sliceDown+1;

            std::cout << "Slice Down " << sliceDown << std::endl;
            std::cout << "Slice Up " << sliceDown +1 << std::endl;

            std::cout << spacing [2] << std::endl;
            std::cout << baseSpacing[2] << std::endl;

            double lambda2= j*spacing[2]/(double)baseSpacing[2]-sliceDown;
            double lambda1=1- lambda2;

            std::cout << "lambda1 " << lambda1 << std::endl;
            std::cout << "lambda2 " << lambda2 << std::endl;

            std::cout << "spacing 2 " << spacing[2] << std::endl;
            std::cout << "basespacing 2 " << baseSpacing[2] << std::endl;


            ExtractSliceFilterType::Pointer extractMaskSlice = ExtractSliceFilterType::New();
            extractMaskSlice->SetInput(maskImage);

            ImageType::RegionType region;

            ImageType::SizeType size=maskImage->GetLargestPossibleRegion().GetSize();
            size[2]=0;
            region.SetSize(size);

            ImageType::IndexType index;
            index.Fill(0);
            index[2]=j;
            region.SetIndex(index);

            extractMaskSlice->SetExtractionRegion(region);
            extractMaskSlice->SetDirectionCollapseToSubmatrix();
            extractMaskSlice->Update();


            SliceImageType::Pointer maskSliceImage = extractMaskSlice->GetOutput();
            maskSliceImage->DisconnectPipeline();



            SliceImageType::SpacingType sliceSpacing;
            sliceSpacing[0]=baseSpacing[0];
            sliceSpacing[1]=baseSpacing[1];

            SliceImageType::PointType sliceOrigin;
            sliceOrigin.Fill(0);

            SliceImageType::IndexType sliceIndex;
            sliceIndex.Fill(0);


            SliceImageType::SizeType sliceSize;
            sliceSize[0]=size[0];
            sliceSize[1]=size[1];

            SliceImageType::RegionType region1;
            region1.SetSize(sliceSize);
            region1.SetIndex(sliceIndex);

            SliceImageType::RegionType region2;
            region2.SetSize(sliceSize);
            region2.SetIndex(sliceIndex);

            SliceImageType::RegionType region3;
            region3.SetSize(sliceSize);
            region3.SetIndex(sliceIndex);

            SliceImageType::DirectionType direction;
            direction.SetIdentity();



            SliceImageType::Pointer downImage = SliceImageType::New();
            downImage->SetOrigin(sliceOrigin);
            downImage->SetSpacing(sliceSpacing);
            downImage->SetRegions(region1);
            downImage->SetDirection(direction);

            downImage->Allocate();

            SliceImageType::Pointer upImage = SliceImageType::New();
            upImage->SetOrigin(sliceOrigin);
            upImage->SetSpacing(sliceSpacing);
            upImage->SetRegions(region2);
            upImage->SetDirection(direction);

            upImage->Allocate();

            SliceImageType::Pointer currentImage = SliceImageType::New();
            currentImage->SetOrigin(sliceOrigin);
            currentImage->SetSpacing(sliceSpacing);
            currentImage->SetRegions(region3);
            currentImage->SetDirection(direction);

            currentImage->Allocate();

            maskSliceImage->SetOrigin(sliceOrigin);
            maskSliceImage->SetSpacing(sliceSpacing);
            maskSliceImage->SetDirection(direction);



            itk::ImageRegionIterator<SliceImageType> itDown(downImage,downImage->GetLargestPossibleRegion());
            itk::ImageRegionIterator<SliceImageType> itUp(upImage,upImage->GetLargestPossibleRegion());
            itk::ImageRegionIterator<SliceImageType> itCurrent(currentImage,currentImage->GetLargestPossibleRegion());

            itk::ImageRegionIterator<ImageType> it(image,image->GetLargestPossibleRegion());


            ImageType::IndexType indexIterator;
            indexIterator.Fill(0);
            indexIterator[2]=sliceDown;

            it.SetIndex(indexIterator);
            itDown.GoToBegin();
            itUp.GoToBegin();
            itCurrent.GoToBegin();

            for (;!itDown.IsAtEnd();++it,++itDown)
            {
                itDown.Set(it.Get());
                if (lambda1>0.5)
                {
                    itCurrent.Set(it.Get());
                    ++itCurrent;
                }
            }

            for (;!itUp.IsAtEnd();++it,++itUp)
            {
                itUp.Set(it.Get());
                if (!(lambda1>0.5))
                {
                    itCurrent.Set(it.Get());
                    ++itCurrent;
                }
            }



            for (int iteration=1; iteration<args.nbrIterations+1;++iteration)
            {

                std::cout << "Iteration number " << iteration << std::endl;

                typedef rpi::LCClogDemons< SliceImageType, SliceImageType, double >
                        RegistrationMethod;

                typedef itk::Transform<double, 2, 2>
                        LinearTransformType;


                typedef itk::Transform< double, 2, 2 >
                        TransformType;

         RegistrationMethod::UpdateRule updateRule;

                switch( args.updateRule )
                           {
                case 0:
                               updateRule=RegistrationMethod::UPDATE_LOG_DOMAIN ; break;
                case 1:
                    updateRule= RegistrationMethod::UPDATE_SYMMETRIC_LOG_DOMAIN ;break;
                case 2:
                    updateRule=  RegistrationMethod::UPDATE_SYMMETRIC_LOCAL_LOG_DOMAIN; break;
                default:
                    throw std::runtime_error( "Update rule must fit in the range [0,2]." );}


                // Creation of the registration object
                RegistrationMethod * registration1 = new RegistrationMethod();


                registration1->SetSigmaI(args.sigmaI);

                registration1->SetStationaryVelocityFieldStandardDeviation(args.sigmaVel);
                registration1->SetRegularizationType(args.regularization);
                registration1->SetUpdateRule(updateRule);
                if (!args.maskImage.empty())
                {
                    registration1->SetMaskImage(maskSliceImage);
                    registration1->UseMask(true);
                }
                 registration1->SetMovingImage(currentImage);
                 registration1->SetFixedImage(downImage);
                 registration1->SetVerbosity(false);
                 registration1->StartRegistration();

                 TransformType * transformp1= registration1->GetTransformation();
                 SliceFieldTransformType * transform1 = dynamic_cast<SliceFieldTransformType *>(transformp1);

                 svf1 = transform1->GetParametersAsVectorField();


                 RegistrationMethod * registration2 = new RegistrationMethod();
                 registration2->SetSigmaI(args.sigmaI);
                 registration2->SetStationaryVelocityFieldStandardDeviation(args.sigmaVel);
                 registration2->SetRegularizationType(args.regularization);
                 registration2->SetUpdateRule(updateRule);

                 registration2->SetMovingImage(currentImage);
                 registration2->SetFixedImage(upImage);
                 registration2->SetVerbosity(false);

                 if (!args.maskImage.empty())
                 {
                     registration2->SetMaskImage(maskSliceImage);
                     registration2->UseMask(true);
                 }

                 registration2->StartRegistration();

                 TransformType * transformp2= registration2->GetTransformation();
                 SliceFieldTransformType * transform2 = dynamic_cast<SliceFieldTransformType *>(transformp2);

                 svf2 = transform2->GetParametersAsVectorField();


                 SliceFieldType::Pointer svfLambda= SliceFieldType::New();
                 svfLambda->CopyInformation(svf1);
                 svfLambda->SetRegions(svf1->GetLargestPossibleRegion());
                 svfLambda->Allocate();


                 itk::ImageRegionIterator<const SliceFieldType> it1(svf1, svf1->GetLargestPossibleRegion());
                 itk::ImageRegionIterator<const SliceFieldType> it2(svf2, svf2->GetLargestPossibleRegion());
                 itk::ImageRegionIterator<SliceFieldType> itl(svfLambda, svfLambda->GetLargestPossibleRegion());

                 std::cout << " Computing lambda weighted sum of svf" << std::endl;

                 for (it1.GoToBegin(),it2.GoToBegin(),itl.GoToBegin(); !it1.IsAtEnd(); ++it1,++it2,++itl)
                 {
                     SliceVectorType newVector=lambda1*it1.Get()+lambda2*it2.Get();
                     itl.Set(newVector);

                 }
           std::cout << " Computing exponential of svf" << std::endl;

                         ExponentialFilterPointer filter= ExponentialFilterType::New();
                         filter->SetInput(svfLambda);
                         //filter->Update();


                         std::cout << " Warping Image" << std::endl;
                         typename WarperType::Pointer warper = WarperType::New();

                         warper->SetInput( currentImage );
                         warper->SetOutputSpacing( currentImage->GetSpacing() );
                         warper->SetOutputOrigin( currentImage->GetOrigin() );
                         warper->SetOutputDirection( currentImage->GetDirection() );
                         warper->SetDeformationField( filter->GetOutput() );

                         warper->Update();

                         currentImage=warper->GetOutput();
                         currentImage->DisconnectPipeline();

                }

            ExponentialFilterPointer filter1= ExponentialFilterType::New();
            filter1->SetInput(svf1);
            filter1->SetComputeInverse(true);
            filter1->Update();

            itk::ImageRegionIterator<const SliceFieldType> it1(svf1, svf1->GetLargestPossibleRegion());
            itk::ImageRegionIterator<const SliceFieldType> it2(svf2, svf2->GetLargestPossibleRegion());

            double norm1=0;
            double norm2=0;

            for (it1.GoToBegin(),it2.GoToBegin(); !it1.IsAtEnd(); ++it1,++it2)
            {
                norm1+=it1.Get().GetNorm();
                norm2+=it2.Get().GetNorm();
            }

            std::cout << "norm 1 " << norm1 << std::endl;
            std::cout << "norm 2 " << norm2 << std::endl;


            typename WarperType::Pointer warper1 = WarperType::New();

            warper1->SetInput( downImage );
            warper1->SetOutputSpacing( downImage->GetSpacing() );
            warper1->SetOutputOrigin( downImage->GetOrigin() );
            warper1->SetOutputDirection( downImage->GetDirection() );
            warper1->SetDeformationField( filter1->GetOutput() );

            warper1->Update();

            SliceImageType::Pointer downImage2=warper1->GetOutput();
            downImage2->DisconnectPipeline();

            ExponentialFilterPointer filter2= ExponentialFilterType::New();
            filter2->SetInput(svf2);
            filter2->SetComputeInverse(true);
            filter2->Update();


            typename WarperType::Pointer warper2 = WarperType::New();

            warper2->SetInput( upImage );
            warper2->SetOutputSpacing( upImage->GetSpacing() );
            warper2->SetOutputOrigin( upImage->GetOrigin() );
            warper2->SetOutputDirection( upImage->GetDirection() );
            warper2->SetDeformationField( filter2->GetOutput() );

            warper2->Update();

            SliceImageType::Pointer upImage2=warper2->GetOutput();
            upImage2->DisconnectPipeline();



            itk::ImageRegionIterator<ImageType> itResampled(resampledImage,resampledImage->GetLargestPossibleRegion());

            ImageType::IndexType indexResampled;
            indexResampled.Fill(0);
            indexResampled[2]=j;

            itResampled.SetIndex(indexResampled);

            itk::ImageFileWriter<SliceImageType>::Pointer writerSliceDown= itk::ImageFileWriter<SliceImageType>::New();

            std::ostringstream nameDown;
            nameDown << args.outputFolder;
            nameDown<<"image_down_slice";
            nameDown<<j<<".mha";
            writerSliceDown->SetFileName(nameDown.str().c_str());
            writerSliceDown->SetInput(downImage2);
            //writerSliceDown->Update();

            itk::ImageFileWriter<SliceImageType>::Pointer writerSliceUp= itk::ImageFileWriter<SliceImageType>::New();

            std::ostringstream nameUp;
            nameUp << args.outputFolder;
            nameUp<<"image_up_slice";
            nameUp<<j<<".mha";
            writerSliceUp->SetFileName(nameUp.str().c_str());
            writerSliceUp->SetInput(upImage2);
            //writerSliceUp->Update();



            itk::ImageRegionIterator<SliceImageType> itDown2(downImage2,downImage2->GetLargestPossibleRegion());
            itk::ImageRegionIterator<SliceImageType> itUp2(upImage2,upImage2->GetLargestPossibleRegion());


            for (itUp2.GoToBegin(),itDown2.GoToBegin();!itUp2.IsAtEnd();++itUp2,++itResampled,++itDown2)
            {
                double pixelValue=lambda1*itDown2.Get()+lambda2*itUp2.Get();
                itResampled.Set(pixelValue);
            }



        }






       ImageFileWriter::Pointer writerImage = ImageFileWriter::New();
       //writerImage->SetFileName(outputNamesGenerator->GetFileNames()[fni]);
       writerImage->SetFileName(args.outputFolder);

       writerImage->SetInput(resampledImage);
       writerImage->Update();
    }
}




