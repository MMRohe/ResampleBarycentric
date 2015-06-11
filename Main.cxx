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

#include "rpiLCClogDemons.hxx"


struct arguments
{

    int          nbrPatient;                        /* -P option */
    int          nbrFrame;                    /* -F option */
    int nbrIterations;                       /* -t option */
    int updateRule;                         /* -r option */
    int regularization;                     /* -R option */
    double sigmaI;                          /* -S option */
    double sigmaVel;                        /* -d option */
    std::string outputFolder;                   /* -o option */
    std::string inputFolder;                   /* -i option */
    std::string inputLambda;                   /* -l option */


    friend std::ostream& operator<< (std::ostream& o, const arguments& args)
    {

        return o
        <<"  Arguments structure:"<<std::endl
        <<"  Patient: "<<args.nbrPatient<<std::endl
        <<"  Frame: "<<args.nbrFrame<<std::endl
        <<"  Number Iterations: "<<args.nbrIterations<<std::endl
        <<" Input Lambda:" << args.inputLambda << std::endl
        <<" Update Rule:" << args.updateRule<< std::endl
        <<" Regularization:" << args.regularization<< std::endl
        <<" Sigma I:" << args.sigmaI<< std::endl
        <<" Sigma Vel:" << args.sigmaVel<< std::endl
        <<" Input Folder:" << args.inputFolder << std::endl
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
    command.SetOption("Patient","P",true,"Patient");
    command.SetOptionLongTag("Patient","input-patient");
    command.AddOptionField("Patient","intval",MetaCommand::STRING,true);

    command.SetOption("Frame","F",true,"Frame");
    command.SetOptionLongTag("Frame","frame");
    command.AddOptionField("Frame","intval",MetaCommand::INT,true);

    command.SetOption("Iteration","t",false,"Iteration");
    command.SetOptionLongTag("Iteration","nbr-Iteration");
    command.AddOptionField("Iteration","intval",MetaCommand::INT,false,"5");

    command.SetOption("OutputFolder","o",false,"Ouput folder");
    command.SetOptionLongTag("OutputFolder","output-folder");
    command.AddOptionField("OutputFolder","filename",MetaCommand::STRING,false,"");

    command.SetOption("InputFolder","i",false,"Input folder");
    command.SetOptionLongTag("InputFolder","input-folder");
    command.AddOptionField("InputFolder","filename",MetaCommand::STRING,false,"");


    command.SetOption("InputLambda","l",false,"Input lambda");
    command.SetOptionLongTag("InputLambda","input-lambda");
    command.AddOptionField("InputLambda","filename",MetaCommand::STRING,false,"");

    command.SetOption("UpdateRule","r",false,"Update rule:  0 : SSD similarity - exp(v)<-exp(v)oexp(u) (log-domain),1 : SSD similarity - exp(v)<-symmetrized(exp(v)oexp(u)) (symmetriclog-domain), 2 : LCC similarity (default 2).");
    command.SetOptionLongTag("UpdateRule","update-rule");
    command.AddOptionField("UpdateRule","intval",MetaCommand::INT,false,"1");

    command.SetOption("regularization","R",false,"Type of Regularization: 0 = Gaussian convolution (classical Demons), 1 = Harmonic + Bending Energy (default 1)");
    command.SetOptionLongTag("regularization","type-reg");
    command.AddOptionField("regularization","intval",MetaCommand::INT,false,"0");


    command.SetOption("SigmaI","S",false,"Trade-off between similarity and regularization (sigma_i ^2): 0 (sharper but unregular deformations) < sigma_i <= 1 (smootherdeformations but weaker correspondencies). Default: 0.15");
    command.SetOptionLongTag("SigmaI","sigma-I");
    command.AddOptionField("SigmaI","floatval",MetaCommand::FLOAT,false,"0.5");

    command.SetOption("SigmaVel","d",false,"Standard deviation of the Gaussian smoothing of the stationary velocity field (world units). Setting it below 0.1 means no smoothing will be performed (default 1.5).");
    command.SetOptionLongTag("SigmaVel","sigma-vel");
    command.AddOptionField("SigmaVel","floatval",MetaCommand::FLOAT,false,"4");





    // Actually parse the command line
    if (!command.Parse(argc,argv))
    {
        exit( EXIT_FAILURE );
    }


    // Store the parsed information into a struct

    args.nbrPatient=command.GetValueAsInt("Patient","intval");
    args.nbrFrame=command.GetValueAsInt("Frame","intval");
    args.outputFolder = command.GetValueAsString("OutputFolder","filename");
    args.inputFolder = command.GetValueAsString("InputFolder","filename");
    args.inputLambda = command.GetValueAsString("InputLambda","filename");
    args.nbrIterations = command.GetValueAsInt("Iteration","intval");
    args.updateRule= command.GetValueAsInt("UpdateRule","intval");
    args.regularization= command.GetValueAsInt("regularization","intval");
    args.sigmaI= command.GetValueAsFloat("SigmaI","floatval");
    args.sigmaVel= command.GetValueAsFloat("SigmaVel","floatval");

}




int main( int argc, char *argv[] )
{
    typedef float TScalar;
    const int Dimension = 3;
    typedef itk::Vector<double,Dimension>                    VectorType;
    typedef itk::Image<VectorType,Dimension>                FieldType;
    typedef itk::ImageFileReader<FieldType>                FieldTypeReader;
    typedef itk::ImageFileWriter<FieldType>                FieldTypeWriter;


    typedef itk::Image<double,Dimension>                   ImageType;
    typedef ImageType::Pointer                              ImagePointer;


    typedef itk::Image<double,Dimension>                   ImageMaskType;
        typedef ImageMaskType::Pointer                              ImageMaskPointer;

    typedef itk::ImageFileReader<ImageMaskType>       ImageMaskReaderType;
    typedef ImageMaskReaderType::Pointer                        ImageMaskReaderPointer;

    typedef itk::Image<float,Dimension>                   ImageLabelType;
    typedef ImageLabelType::Pointer                              ImageLabelPointer;

    typedef itk::ImageFileWriter<ImageType>                   ImageFileWriter;
        typedef ImageFileWriter::Pointer                              ImageFileWriterPointer;

    typedef itk::ImageFileWriter<ImageLabelType>                   ImageLabelWriterType;
        typedef ImageLabelWriterType::Pointer                              ImageLabelWriterPointer;

    typedef itk::ImageFileReader<ImageType>       ImageReaderType;
    typedef ImageReaderType::Pointer                        ImageReaderPointer;

    typedef itk::ImageFileReader<ImageLabelType>       ImageReaderLabelType;
    typedef ImageReaderLabelType::Pointer                        ImageReaderLabelPointer;



    typedef itk::StationaryVelocityFieldTransform<double, 3>  FieldTransformType;

    typedef itk::ExponentialDeformationFieldImageFilter<FieldType,FieldType> ExponentialFilterType;
    typedef ExponentialFilterType::Pointer                                   ExponentialFilterPointer;

    typedef itk::WarpImageFilter< ImageType, ImageType, FieldType >  WarperType;

    struct arguments args;
    parseOpts (argc, argv, args);

    /// Vector type
    typedef vnl_vector<TScalar> VNLVectorType;
    /// Matrix type.
    typedef vnl_matrix<TScalar> VNLMatrixType;
    /// List of matrices type.
    typedef std::vector<VNLMatrixType> VNLMatrixList;

    std::cout << "Starting Lambda registration with the following arguments:" << std::endl;
    std::cout<<args<<std::endl<<std::endl;


     std::ostringstream baseFolderPatiento;
     baseFolderPatiento << args.inputFolder;
     baseFolderPatiento << "Patient_";
     baseFolderPatiento<< args.nbrPatient;

    std::string baseFolderPatient=baseFolderPatiento.str();

    ImageReaderPointer reader = ImageReaderType::New();

    std::ostringstream imageFilename;
    imageFilename<<baseFolderPatient;
    imageFilename<<"/CropImage/Image_frame_";
    imageFilename<<args.nbrFrame;
    imageFilename<<".mha";

    reader->SetFileName(imageFilename.str().c_str());
   reader->Update();

   ImagePointer baseImage= ImageType::New();
   baseImage= reader->GetOutput();


   std::ostringstream lambdaFilename;
   lambdaFilename << args.inputLambda;

   vnl_matrix< float > Lambda;

   std::ifstream fileLambda(lambdaFilename.str().c_str());
   if (fileLambda.is_open()) { /* ok, proceed with output */
       Lambda.read_ascii(fileLambda);
       fileLambda.close();
   }

   std::cout << Lambda << std::endl;

   std::cout << "Number of Rows: " << Lambda.rows() << std::endl;
   std::cout << "Number of Cols: " <<Lambda.cols() << std::endl;

   double lambda1 = Lambda(args.nbrFrame,0);
   double lambda2 = Lambda(args.nbrFrame,1);
   double lambda3 = Lambda(args.nbrFrame,2);

   int ref1 = (int) Lambda(0,0);
   int ref2 = (int) Lambda(0,1);
   int ref3 = (int) Lambda(0,2);

   std::cout << "1st reference is " << ref1 << std::endl;
   std::cout << "2nd reference is " << ref2 << std::endl;
   std::cout << "3rd reference is " << ref3 << std::endl;

   int baseRef;

   if (lambda1>lambda2)
       if (lambda1>lambda3)
           baseRef=ref1;
        else
           baseRef=ref3;
   else
       if (lambda2>lambda3)
           baseRef=ref2;
   else
           baseRef=ref3;

    std::cout << "Base reference is " << baseRef << std::endl;

   std::ostringstream currentFilename;
   currentFilename<<baseFolderPatient;
   currentFilename<<"/CropImage/Image_frame_";
   currentFilename<<baseRef;
   currentFilename<<".mha";


   ImageReaderPointer readerc = ImageReaderType::New();

   readerc->SetFileName(currentFilename.str().c_str());
   readerc->Update();

   ImagePointer currentImage= ImageType::New();
   currentImage= readerc->GetOutput();


   std::cout << "Reading Ref Images " << std::endl;

   std::ostringstream ref1Filename;
   ref1Filename<<baseFolderPatient;
   ref1Filename<<"/CropImage/Image_frame_";
   ref1Filename<<ref1;
   ref1Filename<<".mha";

   std::cout << "First reference is " << ref1Filename.str() << std::endl;

   ImageReaderPointer reader1 = ImageReaderType::New();

   reader1->SetFileName(ref1Filename.str().c_str());
   reader1->Update();

   ImagePointer ImageRef1= ImageType::New();
   ImageRef1= reader1->GetOutput();


   std::ostringstream ref2Filename;
   ref2Filename<<baseFolderPatient;
   ref2Filename<<"/CropImage/Image_frame_";
   ref2Filename<<ref2;
   ref2Filename<<".mha";
   std::cout << "2nd reference is " << ref2Filename.str() << std::endl;

   ImageReaderPointer reader2 = ImageReaderType::New();

   reader2->SetFileName(ref2Filename.str().c_str());
   reader2->Update();


   ImagePointer ImageRef2= ImageType::New();
   ImageRef2= reader2->GetOutput();


   std::ostringstream ref3Filename;
   ref3Filename<<baseFolderPatient;
   ref3Filename<<"/CropImage/Image_frame_";
   ref3Filename<<ref3;
   ref3Filename<<".mha";
   std::cout << "3rd reference is " << ref3Filename.str() << std::endl;

   ImageReaderPointer reader3 = ImageReaderType::New();

   reader3->SetFileName(ref3Filename.str().c_str());
   reader3->Update();


   ImagePointer ImageRef3= ImageType::New();
   ImageRef3= reader3->GetOutput();





   std::ostringstream imageFilename0;
   imageFilename0<<args.outputFolder;
   imageFilename0<<"ImageOutput_I0";
   imageFilename0<<"_F";
   imageFilename0<<args.nbrFrame;
   imageFilename0<<".mha";


   ImageFileWriter::Pointer writerImage0 = ImageFileWriter::New();
   writerImage0->SetFileName(imageFilename0.str().c_str());
   writerImage0->SetInput(currentImage);
   writerImage0->Update();




   const FieldType * svf3= FieldType::New();
   const FieldType * svf2= FieldType::New();
   const FieldType * svf1= FieldType::New();





   for (int iteration=1; iteration<args.nbrIterations+1;++iteration)
   {

       typedef rpi::LCClogDemons< ImageType, ImageType, double >
               RegistrationMethod;

       typedef itk::Transform<double, 3, 3>
               LinearTransformType;


       typedef itk::Transform< double, 3, 3 >
               TransformType;

RegistrationMethod::UpdateRule updateRule;

       switch( args.updateRule )
                  {
       case 0:
                      updateRule=RegistrationMethod::UPDATE_LOG_DOMAIN ; break;
       case 1:
           updateRule=RegistrationMethod::UPDATE_LOG_DOMAIN ;break;
       case 2:
           updateRule= RegistrationMethod::UPDATE_SYMMETRIC_LOG_DOMAIN ;break;
       case 3:
           updateRule=  RegistrationMethod::UPDATE_SYMMETRIC_LOCAL_LOG_DOMAIN; break;
       default:
           throw std::runtime_error( "Update rule must fit in the range [0,2]." );}


       // Creation of the registration object
       RegistrationMethod * registration1 = new RegistrationMethod();


       registration1->SetSigmaI(args.sigmaI);

       registration1->SetStationaryVelocityFieldStandardDeviation(args.sigmaVel);
       registration1->SetRegularizationType(args.regularization);
       registration1->SetUpdateRule(updateRule);


        registration1->SetMovingImage(currentImage);
        registration1->SetFixedImage(ImageRef1);
        registration1->StartRegistration();


        TransformType * transformp1= registration1->GetTransformation();
        FieldTransformType * transform1 = dynamic_cast<FieldTransformType *>(transformp1);

        svf1 = transform1->GetParametersAsVectorField();


        RegistrationMethod * registration2 = new RegistrationMethod();
        registration2->SetSigmaI(args.sigmaI);
        registration2->SetStationaryVelocityFieldStandardDeviation(args.sigmaVel);
        registration2->SetRegularizationType(args.regularization);
        registration2->SetUpdateRule(updateRule);

        registration2->SetMovingImage(currentImage);
        registration2->SetFixedImage(ImageRef2);
        registration2->StartRegistration();

        TransformType * transformp2= registration2->GetTransformation();
        FieldTransformType * transform2 = dynamic_cast<FieldTransformType *>(transformp2);


        svf2 = transform2->GetParametersAsVectorField();

        RegistrationMethod * registration3 = new RegistrationMethod();

        registration3->SetSigmaI(args.sigmaI);
        registration3->SetStationaryVelocityFieldStandardDeviation(args.sigmaVel);
        registration3->SetRegularizationType(args.regularization);
        registration3->SetUpdateRule(updateRule);
        registration3->SetMovingImage(currentImage);
        registration3->SetFixedImage(ImageRef3);
        registration3->StartRegistration();

        TransformType * transformp3= registration3->GetTransformation();
        FieldTransformType * transform3 = dynamic_cast<FieldTransformType *>(transformp3);

        svf3 = transform3->GetParametersAsVectorField();




        FieldType::Pointer svfLambda= FieldType::New();
        svfLambda->CopyInformation(svf1);
        svfLambda->SetRegions(svf1->GetLargestPossibleRegion());
        svfLambda->Allocate();


        itk::ImageRegionIterator<const FieldType> it1(svf1, svf1->GetLargestPossibleRegion());
        itk::ImageRegionIterator<const FieldType> it2(svf2, svf2->GetLargestPossibleRegion());
        itk::ImageRegionIterator<const FieldType> it3(svf3, svf3->GetLargestPossibleRegion());
        itk::ImageRegionIterator<FieldType> itl(svfLambda, svfLambda->GetLargestPossibleRegion());

        std::cout << " Computing lambda weighted sum of svf" << std::endl;

        for (it1.GoToBegin(),it2.GoToBegin(),it3.GoToBegin(),itl.GoToBegin(); !it1.IsAtEnd(); ++it1,++it2,++it3,++itl)
        {
            VectorType newVector=lambda1*it1.Get()+lambda2*it2.Get()+lambda3*it3.Get();
            itl.Set(newVector);

        }

        std::cout << " Computing exponential of svf" << std::endl;

                ExponentialFilterPointer filter= ExponentialFilterType::New();
                filter->SetInput(svfLambda);
                filter->SetMaximumNumberOfIterations(100);
                filter->Update();


                std::cout << " Warping Image" << std::endl;
                typename WarperType::Pointer warper = WarperType::New();

                warper->SetInput( currentImage );
                warper->SetOutputSpacing( currentImage->GetSpacing() );
                warper->SetOutputOrigin( currentImage->GetOrigin() );
                warper->SetOutputDirection( currentImage->GetDirection() );
                warper->SetDeformationField( filter->GetOutput() );
                warper->Update();

                currentImage=warper->GetOutput();

        std::ostringstream svfFilename;
        svfFilename<<args.outputFolder;
        svfFilename<<"SVFOutput_I";
        svfFilename<<iteration;
        svfFilename<<"_F";
        svfFilename<<args.nbrFrame;
        svfFilename<<".mha";

        std::ostringstream imageFilename;
        imageFilename<<args.outputFolder;
        imageFilename<<"ImageOutput_I";
        imageFilename<<iteration;
        imageFilename<<"_F";
        imageFilename<<args.nbrFrame;
        imageFilename<<".mha";




        FieldTypeWriter::Pointer writerSvf = FieldTypeWriter::New();
        writerSvf->SetFileName(svfFilename.str().c_str());
        writerSvf->SetInput(svfLambda);
        //writerSvf->Update();

        ImageFileWriter::Pointer writerImage = ImageFileWriter::New();
        writerImage->SetFileName(imageFilename.str().c_str());
        writerImage->SetInput(currentImage);
        //writerImage->Update();

        if (iteration==args.nbrIterations)
        {
        std::ostringstream svf1Final;
        svf1Final<<args.outputFolder;
        svf1Final<<"SVFOutput_Finalv1";
        svf1Final<<"_F";
        svf1Final<<args.nbrFrame;
        svf1Final<<".mha";

        std::ostringstream svf2Final;
        svf2Final<<args.outputFolder;
        svf2Final<<"SVFOutput_Finalv2";
        svf2Final<<"_F";
        svf2Final<<args.nbrFrame;
        svf2Final<<".mha";

        std::ostringstream svf3Final;
        svf3Final<<args.outputFolder;
        svf3Final<<"SVFOutput_Finalv3";
        svf3Final<<"_F";
        svf3Final<<args.nbrFrame;
        svf3Final<<".mha";





        FieldTypeWriter::Pointer writerSvf = FieldTypeWriter::New();
        writerSvf->SetFileName(svf1Final.str().c_str());
        writerSvf->SetInput(svf1);
        writerSvf->Update();

        writerSvf->SetFileName(svf2Final.str().c_str());
        writerSvf->SetInput(svf2);
        writerSvf->Update();

        writerSvf->SetFileName(svf3Final.str().c_str());
        writerSvf->SetInput(svf3);
        writerSvf->Update();

        }

   }

   ImageRef1->DisconnectPipeline();
   ImageRef2->DisconnectPipeline();
   ImageRef3->DisconnectPipeline();


   ExponentialFilterPointer filter1= ExponentialFilterType::New();
   filter1->SetInput(svf1);
   filter1->ComputeInverseOn();
   filter1->SetMaximumNumberOfIterations(100);
   filter1->Update();

   typename WarperType::Pointer warper1 = WarperType::New();

   warper1->SetInput( ImageRef1 );
   warper1->SetOutputSpacing( ImageRef1->GetSpacing() );
   warper1->SetOutputOrigin( ImageRef1->GetOrigin() );
   warper1->SetOutputDirection( ImageRef1->GetDirection() );
   warper1->SetDeformationField( filter1->GetOutput() );
   warper1->Update();


   ExponentialFilterPointer filter2= ExponentialFilterType::New();
   filter2->SetInput(svf2);
   filter2->ComputeInverseOn();
   filter2->SetMaximumNumberOfIterations(100);
   filter2->Update();


   typename WarperType::Pointer warper2 = WarperType::New();

   warper2->SetInput( ImageRef2 );
   warper2->SetOutputSpacing( ImageRef2->GetSpacing() );
   warper2->SetOutputOrigin( ImageRef2->GetOrigin() );
   warper2->SetOutputDirection( ImageRef2->GetDirection() );
   warper2->SetDeformationField( filter2->GetOutput() );
   warper2->Update();

   ExponentialFilterPointer filter3= ExponentialFilterType::New();
   filter3->SetInput(svf3);
   filter3->ComputeInverseOn();
   filter3->SetMaximumNumberOfIterations(100);
   filter3->Update();


   typename WarperType::Pointer warper3 = WarperType::New();

   warper3->SetInput( ImageRef3 );
   warper3->SetOutputSpacing( ImageRef3->GetSpacing() );
   warper3->SetOutputOrigin( ImageRef3->GetOrigin() );
   warper3->SetOutputDirection( ImageRef3->GetDirection() );
   warper3->SetDeformationField( filter3->GetOutput() );
   warper3->Update();

   ImageType::Pointer imageFinal1=ImageType::New();
   ImageType::Pointer imageFinal2=ImageType::New();
   ImageType::Pointer imageFinal3=ImageType::New();

   imageFinal1=warper1->GetOutput();
   imageFinal2=warper2->GetOutput();
   imageFinal3=warper3->GetOutput();

   ImageType::Pointer imageFinalLambda=ImageType::New();
   imageFinalLambda->CopyInformation(imageFinal1);
   imageFinalLambda->SetRegions(imageFinal1->GetLargestPossibleRegion());
   imageFinalLambda->Allocate();



   itk::ImageRegionIterator<ImageType> it1(imageFinal1, imageFinal1->GetLargestPossibleRegion());
   itk::ImageRegionIterator<ImageType> it2(imageFinal2, imageFinal2->GetLargestPossibleRegion());
   itk::ImageRegionIterator<ImageType> it3(imageFinal3, imageFinal3->GetLargestPossibleRegion());
   itk::ImageRegionIterator<ImageType> itl(imageFinalLambda, imageFinalLambda->GetLargestPossibleRegion());

   std::cout << " Writing image as weighted sum of images" << std::endl;

   for (it1.GoToBegin(),it2.GoToBegin(),it3.GoToBegin(),itl.GoToBegin(); !it1.IsAtEnd(); ++it1,++it2,++it3,++itl)
   {
       ImageType::PixelType newPoint=lambda1*it1.Get()+lambda2*it2.Get()+lambda3*it3.Get();
       itl.Set(newPoint);

   }

   std::ostringstream imageFilenameLambda;
   imageFilenameLambda<<args.outputFolder;
   imageFilenameLambda<<"ImageOutput_Final";
   imageFilenameLambda<<"_F";
   imageFilenameLambda<<args.nbrFrame;
   imageFilenameLambda<<".mha";


   ImageFileWriter::Pointer writerImage = ImageFileWriter::New();
   writerImage->SetFileName(imageFilenameLambda.str().c_str());
   writerImage->SetInput(imageFinalLambda);
   writerImage->Update();


}

