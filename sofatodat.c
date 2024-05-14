#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <netcdf.h>  //netcdf Library used to access data in sofa file
#include <math.h>
#include <fftw3.h>   //Fastest fourier transform in the west library used to interpolate data


#define MAXSTRING 100
#define SQUARE(X) (X)*(X)

char c, filename[MAXSTRING], filenamel[MAXSTRING], filenamer[MAXSTRING],ListenerDimName[MAXSTRING], SourceDimName[MAXSTRING],EmitterDimName[MAXSTRING],ListenerPositionDim[MAXSTRING], SourcePositionDim[MAXSTRING], EmitterPositionDim[MAXSTRING], SourcePositionType[MAXSTRING], SourcePositionUnits[MAXSTRING], SRl[MAXSTRING], SRr[MAXSTRING];
int i = 0;

/* definitions */
/* from mit */
#define minelev (-40)
#define elevincrement (10)
#define truncBits (128)

/* number of measurements per elev: mit data const:read only, static:exists
   for whole process... */
static const int elevationarray[14] =
  {56, 60, 72, 72, 72, 72, 72, 60, 56, 45, 36, 24, 12, 1 };

/*This is the name of the date file to be read */
void stringread(void)
{
  printf("\nSet the location of the SOFA file to be converted:  \n");
  for (i = 0;(c = getchar()) != '\n'; ++i)
  {
    filename[i] = c;
  }
  printf("\n*** file %s chosen for conversion ***\n", filename);


}

/*This is the name of the date file to write */
void stringwrite(double SampleRate)
{
  for (i = 0; i < MAXSTRING; i++)
  {
    filename[i] = 0;
  }
  printf("\nSet the name of the DAT file to be saved  \n");
  printf("\n--this will be appended with *sample rate*l.dat and *sample rate*r.dat for left and right respectively--\n");
  for (i = 0;(c = getchar()) != '\n'; i++)
  {
    filename[i] = c;
  }

  sprintf(SRl,"%i",(int)SampleRate);
  sprintf(SRr,"%i",(int)SampleRate);

  strcpy(filenamel, filename);
  strcpy(filenamer, filename);

  strcat(filenamel,"l");
  strcat(filenamer,"r");

  strcat(SRl,".dat");
  strcat(SRr,".dat");

  strcat(filenamel,SRl);
  strcat(filenamer,SRr);
  printf("\n*** filename %s & %s saving ***\n", filenamel,filenamer);

}

  int main(void)
  {

    FILE *foutl, *foutr, *ValOutTest, *ValDistTest;
    int status = 0 , ncid = 0, ndims = 0, nvars = 0, ngatts = 0, unlimdimid = 0;
    int minimumDelayValue = 0;
    int noofsamplesoutput = 0;
    int samplenumber = 0;

    /*char datatype[MAXSTRING];*/

    int data_ir_id = 0;
    int version_id = 0;
    nc_type vr_type, t_type;
    nc_type l_type, s_type, e_type;
    size_t  vr_len, t_len, EmitterPositionLength[MAXSTRING], SourcePositionTypeLength, SourcePositionUnitsLength, ListenerPositionLength[MAXSTRING], SourcePositionLength[MAXSTRING];

    //Arrays used to define the sizes of Netcdf data structures, used in pulling data from sofa Files
    // Used for 3d array of Samples/Hrtf measurements
    size_t start[3] = {0, 0, 0};
    size_t count[3] = {0, 0, 0};

    //used for listener position Array
    size_t startL[2] = {0, 0};
    size_t countL[2] = {0, 0};

    //Used for source position Array
    size_t startS[2] = {0, 0};
    size_t countS[2] = {0, 0};

    //Used for emitter position Array
    size_t startE[3] = {0, 0, 0};
    size_t countE[3] = {0, 0, 0};

    //Used for Sample Rate Array
    size_t startSR[1] = {0};
    size_t countSR[1] = {0};

    size_t len = 0;

    //Netcdf variables (for Matrixes sizes)
    unsigned long M = 0; //Defines the Number of measurements
    unsigned long R = 0; //Number of recievers (ears sampled)
    unsigned long N = 0; //Number of Samples (per measurement)
    unsigned long C = 0; // Number of Dimensions in Positional data - Most often 3
    unsigned long I = 0;

    //Net cdf meta data to aid in accessing the matrices
    int ListenerPositionID =0, SourcePositionID =0, EmitterPositionID=0;
    int ListenerDimNum =0, SourceDimNum =0, EmitterDimNum =0;
    int ListenerDimID[MAXSTRING], SourceDimID[MAXSTRING], EmitterDimID[MAXSTRING];
    int ListenerAttNum =0, SourceAttNum =0, EmitterAttNum =0;

    //Data IR Meta data tag, and Array to store dim values
    int data_ir_dim = 0;
    int data_dim_vals[10] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

    //For loop variables
    int p = 0, r = 0, s = 0, q = 0, m = 0, j = 0, k = 0, l = 0, jump = 0, jumpb = 0;

    //Arrays and variables to set theoretical values of MIT Kemar data (standard data set used in open source tools)
    double CompareSourcePosition[710*3];
    double TestAngleHi = 0, TestAngleLo =0, SPAngle =0;
    float data_HRTF_Outputl[710*truncBits];
    float data_HRTF_Outputr[710*truncBits];
    float *pdata_HRTF_Outputl = NULL, *pdata_HRTF_Outputr = NULL;
    int elevationMIT =0;
    float angle =0;

    //Pointers set up, for later Malloc and Array definition
    double *data_ir_vals = NULL, *data_ir_valsTrunc = NULL, *EmitterPosition = NULL; //3 deep pointer Necessary for 3D Array DynamicMemory allocation
    double *ListenerPosition = NULL, *SourcePosition = NULL, *SpherePosition = NULL; //2 deep pointer to set up 2d array for Dynamic Allocation
    int *truncAddress = NULL; //Setting up for TruncAddress Memory allocation (1d array)

    int SampleRateID =0;
    double data_SampleRate = 0;

    //Variables for spherical positioning
    double azimuth =0, zenith = 0, radialdistance=0;

    //Variables for cartesian source position
    double x=0, y=0, z=0, xa = 0, ya = 0, za = 0;
    double *distance = NULL;
    double min[4] = {0, 0, 0, 0}; //Used for sorting closest values
    double totaldistance = 0;
    double totalweight = 0;
    double weightA = 0, weightB = 0, weightC = 0, weightD = 0;
    double temp = 0;
    int location[4] = {0,1,2,3};


    double interpolAl[truncBits];
    double interpolAr[truncBits];
    double interpolBl[truncBits];
    double interpolBr[truncBits];
    double interpolCl[truncBits];
    double interpolCr[truncBits];
    double interpolDl[truncBits];
    double interpolDr[truncBits];

    double polarAl[truncBits];
    double polarAr[truncBits];
    double polarBl[truncBits];
    double polarBr[truncBits];
    double polarCl[truncBits];
    double polarCr[truncBits];
    double polarDl[truncBits];
    double polarDr[truncBits];

    double fftSampAl[truncBits];
    double fftSampAr[truncBits];
    double fftSampBl[truncBits];
    double fftSampBr[truncBits];
    double fftSampCl[truncBits];
    double fftSampCr[truncBits];
    double fftSampDl[truncBits];
    double fftSampDr[truncBits];

    double HRTFfreql[truncBits], HRTFfreqr[truncBits];
    //Tag for interpolation, and data transfer checked
    int interpolNeed = 0;
    int ignore = 0;
    int NumberofMeasurementsEqual = 0;
    int CarttoSphere = 0;

    //Set up of FFTW3 "plan" variables
    fftw_plan SampAplanl,SampAplanr,SampBplanl,SampBplanr,SampCplanl,SampCplanr,SampDplanl,SampDplanr;


    /* setup fft plans (see fftw documentation) */
  	SampAplanl = fftw_plan_r2r_1d(truncBits, interpolAl, fftSampAl, FFTW_R2HC, FFTW_MEASURE);
    SampAplanr = fftw_plan_r2r_1d(truncBits, interpolAr, fftSampAr, FFTW_R2HC, FFTW_MEASURE);
  	SampBplanl = fftw_plan_r2r_1d(truncBits, interpolBl, fftSampBl, FFTW_R2HC, FFTW_MEASURE);
    SampBplanr = fftw_plan_r2r_1d(truncBits, interpolBr, fftSampBr, FFTW_R2HC, FFTW_MEASURE);
    SampCplanl = fftw_plan_r2r_1d(truncBits, interpolCl, fftSampCl, FFTW_R2HC, FFTW_MEASURE);
    SampCplanr = fftw_plan_r2r_1d(truncBits, interpolCr, fftSampCr, FFTW_R2HC, FFTW_MEASURE);
    SampDplanl = fftw_plan_r2r_1d(truncBits, interpolDl, fftSampDl, FFTW_R2HC, FFTW_MEASURE);
    SampDplanr = fftw_plan_r2r_1d(truncBits, interpolDr, fftSampDr, FFTW_R2HC, FFTW_MEASURE);

    stringread();

    if (access(filename, F_OK) != -1) //Check that the file exists
    {
      //ValDistTest = fopen("distanceTestdata.txt", "wb");

      //**********************
      //Opening the File given by string input using netcdf
      status = nc_open(filename, NC_NOWRITE, &ncid);
        if (status != NC_NOERR)
        {
          printf("\nERROR -- Could not open\n");
          printf("%s\n", nc_strerror(status));
        }

        //Inquiry to the file for basics - number of dimension, number of variables,
        //number of attributes, and the ID for the unlimited Dimension of the DATA
        status = nc_inq(ncid, &ndims, &nvars, &ngatts, &unlimdimid);
        if (status != NC_NOERR)
        {
          printf("\nERROR -- Inquiry Failed\n");
          printf("%s\n", nc_strerror(status));
        }

        //Test for Number of Dimensions and ID tags
        //printf("\nFile ID = %i\n Number of Dimensions = %i\n Number of Variables = %i\n Number of Attributes = %i\n Number of unlimited dimensions = %i\n", ncid, ndims, nvars, ngatts, unlimdimid);


        //**********************
        //Getting the Variable ID for the Data IR
        status = nc_inq_varid(ncid, "Data.IR", &data_ir_id);
        if (status != NC_NOERR)
        {
          printf("\nERROR -- Variable ID failed\n");
          printf("%s\n", nc_strerror(status));
        }
        //Test for Data IR ID
        //printf("\nData IR variable ID = %i\n",data_ir_id);

        //**********************
        /*Finding Global Dimensional Variable M and its length*/
        status = nc_inq_dimid(ncid, "M", &data_ir_dim);
        if (status != NC_NOERR)
        {
          printf("\nERROR -- Variable dimensions failed\n");
          printf("%s\n", nc_strerror(status));
        }

        status = nc_inq_dimlen(ncid, data_ir_dim, &M);
        if (status != NC_NOERR)
        {
          printf("\nERROR -- Variable dimensions failed\n");
          printf("%s\n", nc_strerror(status));
        }
        //Test For M read Correctly
        //printf("\nM eauals %lu\n",M);


        //**********************
        /*Finding Global Dimensional Variable R and its length*/
        status = nc_inq_dimid(ncid, "R", &data_ir_dim);
        if (status != NC_NOERR)
        {
          printf("\nERROR -- Variable dimensions failed\n");
          printf("%s\n", nc_strerror(status));
        }

        status = nc_inq_dimlen(ncid, data_ir_dim, &R);
        if (status != NC_NOERR)
        {
          printf("\nERROR -- Variable dimensions failed\n");
          printf("%s\n", nc_strerror(status));
        }
        //Test for R read Correctly
        //printf("\nR equals %lu\n",R);


        //**********************
        /*Finding Global Dimensional Variable N and its length*/
        status = nc_inq_dimid(ncid, "N", &data_ir_dim);
        if (status != NC_NOERR)
        {
          printf("\nERROR -- Variable dimensions failed\n");
          printf("%s\n", nc_strerror(status));
        }

        status = nc_inq_dimlen(ncid, data_ir_dim, &N);
        if (status != NC_NOERR)
        {
          printf("\nERROR -- Variable dimensions failed\n");
          printf("%s\n", nc_strerror(status));
        }
        //Test for N Read Correctly
        //printf("\nN equals %lu\n",N);

        //**********************
        /*Finding Global Dimensional Variable C and its length*/
        status = nc_inq_dimid(ncid, "C", &data_ir_dim);
        if (status != NC_NOERR)
        {
          printf("\nERROR -- Variable dimensions failed\n");
          printf("%s\n", nc_strerror(status));
        }

        status = nc_inq_dimlen(ncid, data_ir_dim, &C);
        if (status != NC_NOERR)
        {
          printf("\nERROR -- Variable dimensions failed\n");
          printf("%s\n", nc_strerror(status));
        }
        //Test for C read Correctly
        //printf("\nC equals %lu\n",C);

        //**********************
        /*Finding Global Dimensional Variable I and its length*/
        status = nc_inq_dimid(ncid, "I", &data_ir_dim);
        if (status != NC_NOERR)
        {
          printf("\nERROR -- Variable dimensions failed\n");
          printf("%s\n", nc_strerror(status));
        }

        status = nc_inq_dimlen(ncid, data_ir_dim, &I);
        if (status != NC_NOERR)
        {
          printf("\nERROR -- Variable dimensions failed\n");
          printf("%s\n", nc_strerror(status));
        }
        //Test for i read Correctly
        //printf("\nI equals %lu\n",I);

        status = nc_inq_varid(ncid, "Data.SamplingRate", &SampleRateID);
        if (status != NC_NOERR)
        {
          printf("\nERROR -- Sample Rate failed\n");
          printf("%s\n", nc_strerror(status));
        }
        //Test for Data IR ID
        //printf("\nData SampleRate variable ID = %i\n",SampleRateID);

        countSR[0] = I;

        status = nc_get_vara(ncid, SampleRateID, startSR, countSR, &data_SampleRate);
        if (status != NC_NOERR)
        {
          printf("\nERROR -- Get Variable failed\n");
          printf("%s\n", nc_strerror(status));
        }
        //printf("\nData SampleRate = %f\n", data_SampleRate);


        //**********************
        /*Setting up IR Value Matrix Dynamically, and count systems*/
        data_ir_vals = malloc((M*R*N) * sizeof(double));

        if(data_ir_vals == NULL)
        {
          printf("\n\nFailed to allocate data_ir_vals.\n\n");
        }
        count[0] = M;
        count[1] = R;
        count[2] = N;


        //**********************
        //Reading and Storing Data Values
        //As data_ir_vals is already a pointer - & is not necessary
        status = nc_get_vara(ncid, data_ir_id, start, count, data_ir_vals);
        if (status != NC_NOERR)
        {
          printf("\nERROR -- Get Variable failed\n");
          printf("%s\n", nc_strerror(status));
        }


        //TESTING Data IR VALUES BY PRINT
        /*jump = 0;
        for(i=0;i<M;i++) //Cycles through recievers
        {
            for (j=0;j<R;j++) //Cycles through Measurements
            {
              for(k=0;k<N;k++) //Cycles through Samples
                {
                  printf("\n%lf\n\n", data_ir_vals[k+jump]);
                  printf("\nRead at Position %i %i %i %d\n\n", i,j,k, k+jump);
                }

                if (j < (R-1))
                {
                  jump += N;
                }
            }
        jump += N;
        }
        */


//**********************
//**********************
//Reading Listener Data
//**********************
//**********************


        //**********************
        /*Finding ListenerPostion ID*/
        status = nc_inq_varid(ncid, "ListenerPosition", &ListenerPositionID);
        if (status != NC_NOERR)
        {
          printf("\nERROR -- No ID for ListenerPosition\n");
          printf("%s\n", nc_strerror(status));
        }
        //Test for Listener Position ID
        //printf("ListenerPostion ID = %i\n\n",ListenerPositionID);


        //**********************
        //Finding ListenerPosition attributes
        status = nc_inq_var(ncid, ListenerPositionID, 0, &l_type, &ListenerDimNum, &ListenerDimID[0], &ListenerAttNum);
        if (status != NC_NOERR)
        {
          printf("\nERROR -- cannont read ListenerPositionLength\n");
          printf("%s\n", nc_strerror(status));
        }
        //Test For listener Dimension Number
        //printf("\nListenerDimNum = %i", ListenerDimNum);


        //**********************
        //Reading the ID, Name and Length of the dimensions defining Listener Position
        for (i=0;i<ListenerDimNum;i++)
        {
          status = nc_inq_dim(ncid, ListenerDimID[i], &ListenerDimName[i], &ListenerPositionLength[i]);
          if (status != NC_NOERR)
          {
            printf("\nERROR -- cannont read ListenerPositionLength\n");
            printf("%s\n", nc_strerror(status));
          }
          //Testing Attributes by Printing Values
          //printf("\nListenerDimID = %i", ListenerDimID[i]);
          //printf("\nListenerDimName = %c", ListenerDimName[i]);
          //printf("\nListenerPosition Length = %lu\n\n", ListenerPositionLength[i]);
        }


        //printf("\nListenerAttNum = %i\n\n", ListenerAttNum);

        ListenerPosition = malloc((ListenerPositionLength[0]*ListenerPositionLength[1]) * sizeof(double));

        if(ListenerPosition == NULL)
        {
          perror("\n\nFailed to allocate ListenerPosition.\n\n");
        }

        countL[0] = ListenerPositionLength[0];
        countL[1] = ListenerPositionLength[1];

        //Reading Listener Postion
        status = nc_get_vara(ncid, ListenerPositionID, startL, countL, ListenerPosition);
        if (status != NC_NOERR)
        {
          printf("\nERROR -- cannont read ListenerPosition\n");
          printf("%s\n", nc_strerror(status));
        }

        //Testing Listener Position
        /*jump = 0;
        for (i = 0; i< (ListenerPositionLength[0]); i++)
        {
          for(j = 0; j<(ListenerPositionLength[1]); j++)
          {
            printf("\nListenerPosition = %f\n", ListenerPosition[j +jump]);
          }
          jump += ListenerPositionLength[1];
        }
        */


//**********************
//**********************
//Reading Source Data
//**********************
//**********************


        //**********************
        //Finding SourcePosition ID
        status = nc_inq_varid(ncid, "SourcePosition", &SourcePositionID);
        if (status != NC_NOERR)
        {
          printf("\nERROR -- No ID for SourcePosition\n");
          printf("%s\n", nc_strerror(status));
        }

        //Testing Source Position ID
        //printf("\n\nSourcePosition ID = %i\n\n",SourcePositionID);

        //Finding SourcePosition attributes
        status = nc_inq_var(ncid, SourcePositionID, 0, &s_type, &SourceDimNum, &SourceDimID[0], &SourceAttNum);
        if (status != NC_NOERR)
        {
          printf("\nERROR -- cannont read SourcePositionLength\n");
          printf("%s\n", nc_strerror(status));
        }

        //Test read Source Dimension and Attribute Values
        //printf("\nSourceDimNum = %i", SourceDimNum);
        //printf("\nSourceAttNum = %i\n\n", SourceAttNum);


        //**********************
        //Reading and Printing Source Position dimensions
        for (i=0;i<=SourceDimNum;i++)
        {
          status = nc_inq_dim(ncid, SourceDimID[i], &SourceDimName[i], &SourcePositionLength[i]);
          if (status != NC_NOERR)
          {
            printf("\nERROR -- cannont read SourcePositionLength\n");
            printf("%s\n", nc_strerror(status));
          }

          //Test Read Source Position values
          //printf("\nSourceDimID = %i", SourceDimID[i]);
          //printf("\nSourceDimName = %c", SourceDimName[i]);
          //printf("\nSourcePosition Length = %lu\n\n", SourcePositionLength[i]);
        }


        //**********************
        //setting up variables to read Source position dynamically as a 2d array
        SourcePosition = malloc((SourcePositionLength[0]*SourcePositionLength[1]) * sizeof(double));

        if(SourcePosition == NULL)
        {
          perror("\n\nFailed to allocate SourcePosition.\n\n");
        }

        //**********************
        //setting up Array for Spherical position dynamically as a 2d array
        SpherePosition = malloc((SourcePositionLength[0]*SourcePositionLength[1]) * sizeof(double));

        if(SpherePosition == NULL)
        {
          perror("\n\nFailed to allocate SpherePosition.\n\n");
        }

        countS[0] = SourcePositionLength[0];
        countS[1] = SourcePositionLength[1];

        //**********************
        //Reading Source Postion
        status = nc_get_vara(ncid, SourcePositionID, startS, countS, SourcePosition);
        if (status != NC_NOERR)
        {
          printf("\nERROR -- cannont read SourcePosition\n");
          printf("%s\n", nc_strerror(status));
        }


        //Testing Source Position
        /*jump = 0;
        for (i = 0; i < SourcePositionLength[0]; i++)
        {
          for (j = 0; j < SourcePositionLength[1]; j++)
          {
            fprintf(ValDistTest,"\nSourcePosition = %f\n\n", SourcePosition[j + jump]);
          }
          jump += SourcePositionLength[1];
        }
        */

        //**********************
        //Getting Source Position Type Length
        status = nc_inq_attlen(ncid, SourcePositionID, "Type", &SourcePositionTypeLength);
        if (status != NC_NOERR)
        {
          printf("\nERROR -- Cannot Read Length for SourcePosition:Type\n");
          printf("%s\n", nc_strerror(status));
        }

        //Test read Value
        //printf("\n\nSourcePosition:Type Length = %zu\n\n",SourcePositionTypeLength);


        //**********************
        //Getting Source Position Units length
        status = nc_inq_attlen(ncid, SourcePositionID, "Units", &SourcePositionUnitsLength);
        if (status != NC_NOERR)
        {
          printf("\nERROR -- Cannot Read length for SourcePosition:Units\n");
          printf("%s\n", nc_strerror(status));
        }

        //Test Read Value
        //printf("\n\nSourcePosition:Units ID = %zu\n\n",SourcePositionUnitsLength);


        //**********************
        //Reading Source Positional type
        status = nc_get_att(ncid, SourcePositionID,"Type", &SourcePositionType);
        if (status != NC_NOERR)
        {
          printf("\nERROR -- cannont read SourcePosition:Type\n");
          printf("%s\n", nc_strerror(status));
        }

        //**********************
        //Reading Source Positional Units
        status = nc_get_att(ncid, SourcePositionID,"Units", &SourcePositionUnits);
        if (status != NC_NOERR)
        {
          printf("\nERROR -- cannont read SourcePosition:Units\n");
          printf("%s\n", nc_strerror(status));
        }


        //**********************
        //Testing SourcePosition Units
        printf("\nSourcePosition:Type = ");
        for (len = 0; len < SourcePositionTypeLength; len++)
        {
          printf("%c", SourcePositionType[len]);
        }
        printf("\n");


        //Testing SourcePosition Units
        printf("\nSourcePosition:Units = ");
        for (len = 0; len < SourcePositionUnitsLength; len++)
        {
          printf("%c", SourcePositionUnits[len]);
        }
        printf("\n");


//**********************
//**********************
//Reading Emitter Data
//**********************
//**********************


        //**********************
        //Finding EmitterPostion ID
        status = nc_inq_varid(ncid, "EmitterPosition", &EmitterPositionID);
        if (status != NC_NOERR)
        {
          printf("\nERROR -- No ID for EmitterPosition\n");
          printf("%s\n", nc_strerror(status));
        }
        //printf("\n\nEmitterPosition ID = %i\n\n",EmitterPositionID);

        status = nc_inq_var(ncid, EmitterPositionID, 0, &s_type, &EmitterDimNum, &EmitterDimID[0], &EmitterAttNum);
        if (status != NC_NOERR)
        {
          printf("\nERROR -- cannont read EmitterPosition Variables\n");
          printf("%s\n", nc_strerror(status));
        }

        //Test read Values
        //printf("\nEmitterDimNum = %i", EmitterDimNum);
        //printf("\nEmitterAttNum = %i\n\n", EmitterAttNum);


        //**********************
        //Reading and Printing Emitter Position dimensions
        for (i=0;i<=EmitterDimNum;++i)
        {
          status = nc_inq_dim(ncid, EmitterDimID[i], &EmitterDimName[i], &EmitterPositionLength[i]);
          if (status != NC_NOERR)
          {
            printf("\nERROR -- cannont read EmitterPosition Length\n");
            printf("%s\n", nc_strerror(status));
          }

          /*
          //Test Read values
          printf("\nEmitterDimID = %i", EmitterDimID[i]);
          printf("\nEmitterDimName = %c", EmitterDimName[i]);
          printf("\nEmitterPosition Length = %lu\n\n", EmitterPositionLength[i]);
          */
        }


        //**********************
        //setting up variables to read Emitter position
        EmitterPosition = malloc((EmitterPositionLength[0]*EmitterPositionLength[1]*EmitterPositionLength[2]) * sizeof(double));

        if(EmitterPosition == NULL)
        {
          perror("\n\nFailed to allocate EmitterPosition.\n\n");
        }

        countE[0] = EmitterPositionLength[0];
        countE[1] = EmitterPositionLength[1];
        countE[2] = EmitterPositionLength[2];


        //**********************
        //Reading Emitter Postion
        status = nc_get_vara(ncid, EmitterPositionID, startE, countE, EmitterPosition);
        if (status != NC_NOERR)
        {
          printf("\nERROR -- cannont read EmitterPosition\n");
          printf("%s\n", nc_strerror(status));
        }


        //Testing Emitter Position
        /*jump = 0;
        for (int i = 0; i< (EmitterPositionLength[0]); i++)
        {
          for (int j = 0; j< (EmitterPositionLength[1]); j++)
          {
            for (int k = 0; k< (EmitterPositionLength[2]); k++)
            {
              printf("\nEmitterPosition = %f\n", EmitterPosition[k+jump]);
            }
            if(j < (EmitterPositionLength[1]-1))
            {
              jump += EmitterPositionLength[2];
            }
          }
          jump += EmitterPositionLength[2];
        }
        */


        //**********************
        //Checking if Source Position is spherical or Cartesian
        if (strncmp(SourcePositionType,"spherical", SourcePositionTypeLength) || strncmp(SourcePositionType, "Spherical", SourcePositionTypeLength))
        {
          CarttoSphere = 0;
        }
        else
        {
          CarttoSphere = 1;
        }


        //**********************
        //If Cartesian - they will be processed into Spherical
        if (CarttoSphere == 1)
        {
          jump = 0;
          printf("\nHaving been tested Source Position type is identified as Cartesian co ordinates - these will be corrected to Spherical before further processing\n\n");

          for (l = 0; l <= SourcePositionLength[0]; l++)
          {
            x = SourcePosition[l + jump];
            y = SourcePosition[l + jump + 1];
            z = SourcePosition[l + jump + 2];

            //Testing Cartesian Values being read
            //printf("\nx = %f", x);
            //printf("\ny = %f", y);
            //printf("\nz = %f\n\n", z);

            azimuth = atan(y/x) ; //Radians
            azimuth *= 180/M_PI;
            zenith = atan(z/sqrt((SQUARE(x) + SQUARE(y)))); //Radians
            zenith *= 180/M_PI;
            radialdistance = sqrt(SQUARE(x) * SQUARE(y) * SQUARE(z));

            /*
            //Testing Values being read
            printf("\n\nazimuth = %lu", azimuth);
            printf("\n\nazenith = %lu", zenith);
            printf("\n\nradial distance = %lu", radialdistance);
            */

            //Setting Cartesian Position Matrix
            SpherePosition[l + jump] = azimuth;
            SpherePosition[l + jump + 1] = zenith;
            SpherePosition[l + jump + 2] = radialdistance;

            jump += SourcePositionLength[1];
          }

          for (l=0; l<SourcePositionLength[0]*SourcePositionLength[1]; l++)
          {
            SourcePosition[l] = SpherePosition[l];
          }
        }


//**********************
//**********************
//Interpolation Check
//**********************
//**********************



        //**********************
        //Set Parameters of MIT Kemar
        //**********************
        //Calculate and place in Compare Source Array
        s = 0;
        r = 0;
        for (elevationMIT = minelev; elevationMIT <= 90; elevationMIT += elevincrement)
        {
          angle = 360.00000/(elevationarray[s]);
          for (q = 0; q < elevationarray[s]; q++)
          {
            CompareSourcePosition[r] = angle * q;
            CompareSourcePosition[r+1] = elevationMIT;
            CompareSourcePosition[r+2] = 1.4;
            r+=3;
          }
          s++;
        }

        // Test Source Possition issues by Print - Compare Source Positions
        /*jump = 0;
        for (int i = 0; i< (710); i++)
        {
          for (int j = 0; j< (3); j++)
          {
            printf("\nCompareSourcePosition = %f\n", CompareSourcePosition[j + jump]);
            printf("\nSourcePosition = %f\n", SourcePosition[j + jump]);
          }
          jump += 3;
        }
        */


        //**********************
        //Test Number of Measurements against MIT KEMAR Measurement Numbers
        if (M != 710)
        {
          NumberofMeasurementsEqual = 0;
          interpolNeed = 1;
        }
        else
        {
          NumberofMeasurementsEqual = 1;
        }

        if (NumberofMeasurementsEqual == 1)
        {
          //**********************
          //Test Source position from SOFA against Source Position of MIT Calculated.
          jump = 0;
          for (i = 0; i < (SourcePositionLength[0]); i++)
          {
            for (j = 0; j < (SourcePositionLength[1]); j++)
            {
              if (roundf(SourcePosition[j + jump] * 1000)/1000 != roundf(CompareSourcePosition[j + jump] * 1000)/1000)
              {
                //printf("\n\issues at i = %d, j = %d\n\n", i,j);
                //printf("\n\nNeed to Interpolate data\n\n");
                //printf("\n\n Source Position \n%f \n\n", SourcePosition[i][j]);
                //printf("\n\n Compare Position \n%f \n\n", CompareSourcePosition[i][j]);

                //Flag interpolation as needed
                interpolNeed = 1;
                break;
              }
              else
              {
              //Flag interpolation as not needed
              interpolNeed = 0;
              }
            }
            jump +=SourcePositionLength[1];
            if (interpolNeed == 1)
            {
              break;
            }
          }
        }

        //Inform the User
        if (interpolNeed == 1)
        {
          printf("\n\n Interpolation is Needed\n\n");
        }
        else
        {
          printf("\n\n Interpolation is not Needed\n\n");
        }

        //**********************
        //Check if Interpolation is needed and Not equal to 128/Truncated bit length (128 is Mit Kemar)
        if (interpolNeed == 1 && N != truncBits)
        {
          //**********************
          //Allocate a 1D array - to store values "address" of first non Zero value
          truncAddress = malloc(M * sizeof(int));

          if(truncAddress == NULL)
          {
            printf("\n\nFailed to allocate truncAddress.\n\n");
          }

          printf("\n\n Values require Truncation - array set up for earliest address of non 0 values");

          //**********************
          //Set up 3D array to store truncated impulse responses
          data_ir_valsTrunc = malloc((M*R*truncBits) * sizeof(double));

          if(data_ir_valsTrunc == NULL)
          {
            printf("\n\nFailed to allocate data_ir_valsTrunc.\n\n");
          }
          printf("\n\n Data Array set up to store truncated data set\n\n");

          //**********************
          //Test for two points in a row NonZero
          jump = 0;
          for (i=0;i<M;i++) //Step Through Measurements
          {
            for(j=0;j<R;j++) //Step Through recievers
            {
              for(k=0;k<N;k++) //Step through samples
              {
                //**********************
                //Find two NonZero points in a row after the second value
                if(k > 1 && k < (N - 1))
                {
                  if ((data_ir_vals[k + jump] > 0) && (data_ir_vals[k+ jump + 1] > 0))
                  {
                    //Set Starting Point to before the two NonZero Points
                    truncAddress[i] = (k - 1);
                    break;
                  }
                  if ((data_ir_vals[k + jump] < 0) && (data_ir_vals[k+ jump +1]) < 0)
                  {
                    truncAddress[i] = (k - 1);
                    break;
                  }
                }

                //**********************
                //if two nonzero point happen before second value - set First memory location as the memory allocation to start truncating from
                else if(k <= 1 )
                {
                  if (data_ir_vals[k + jump] > 0 && data_ir_vals[k + jump +1] > 0)
                  {
                    //Set Starting Point to before the two NonZero Points
                    truncAddress[i] = 0;
                    break;
                  }
                  else if (data_ir_vals[k + jump] < 0 && data_ir_vals[k + jump +1] < 0)
                  {
                    //Set Starting Point to before the two NonZero Points
                    truncAddress[i] = 0;
                    break;
                  }
                }

                else if(k == (N - 1))
                {
                  if (data_ir_vals[k + jump] > 0 && data_ir_vals[k + jump -1] > 0)
                  {
                    //Set Starting Point to before the two NonZero Points
                    truncAddress[i] = 0;
                    break;
                  }
                  else if (data_ir_vals[k + jump] < 0 && data_ir_vals[k + jump -1] < 0)
                  {
                    //Set Starting Point to before the two NonZero Points
                    truncAddress[i] = 0;
                    break;
                  }

                  }
                }
                if(j < (R-1))
                {
                  jump += N;
                }
              }
              jump +=N;
            }


        //Test Values by Print
        /*for(i=0; i<M;i++)
        {
          printf("\n\nTruncAddress %i = %i", i, truncAddress[i]);
        }
        */


        printf("\n\n Finding earliest address - common for consistent delay\n\n");

        //**********************
        //Find Minimum Value in TruncAddress
        minimumDelayValue = truncBits;
        for (i=0;i<M;i++)
        {
          if(truncAddress[i] < minimumDelayValue)
          {
            minimumDelayValue = truncAddress[i];
          }
        }

        //printf("\n\nAll working fine before Value Transfer\n\n");
        jump = 0;
        jumpb = 0;
        for(i=0;i<M;i++) //Cycles through Measurements
        {
            for (j=0;j<R;j++) //Cycles through Recievers
            {
              for(k=minimumDelayValue;k<(truncBits+minimumDelayValue);k++) //Cycles through truncated Length samples
                {
                  data_ir_valsTrunc[(k-minimumDelayValue)+jump] = data_ir_vals[k+jumpb];
                }

                if (j < (R-1))
                {
                  jumpb += N;
                  jump += truncBits;
                }
            }
        jumpb+= N;
        jump += truncBits;
        }

        //Add fade out to last ten samples
        jump = 0;
        for(i=0;i<M;i++) //Cycles through recievers
        {
            for (j=0;j<R;j++) //Cycles through Measurements
            {
              for(k=0;k<truncBits;k++) //Cycles through Samples
              {
                if(k >= (truncBits - 10) && k < truncBits)
                {
                  data_ir_valsTrunc[k+jump] *= ((truncBits - (k+1))/10);
                }

              }
              if (j < (R-1))
              {
                jump += truncBits;
              }
            }
            jump += truncBits;
        }

      //printf("\n\n Fade out added\n\n");
      }//end of if statement when interpol is needed and number of samples is not equal

      if (interpolNeed == 0 && N != truncBits)
      {
        //**********************
        //Allocate a 1D array - to store values "address" of first non Zero value
        truncAddress = malloc(M * sizeof(int));

        if(truncAddress == NULL)
        {
          printf("\n\nFailed to allocate truncAddress.\n\n");
        }

        printf("\n\n Values require Truncation - array set up for earliest address of non 0 values");

        //**********************
        //Set up 3D array to store truncated impulse responses
        data_ir_valsTrunc = malloc((M*R*truncBits) * sizeof(double));

        if(data_ir_valsTrunc == NULL)
        {
          printf("\n\nFailed to allocate data_ir_valsTrunc.\n\n");
        }
        printf("\n\n Data Array set up to store truncated data set\n\n");

        //**********************
        //Test for two points in a row NonZero
        jump = 0;
        for (i=0;i<M;i++) //Step Through Recievers
        {
          for(j=0;j<R;j++) //Step Through Measurements
          {
            for(k=0;k<N;k++) //Step through samples
            {
              //**********************
              //Find two NonZero points in a row after the second value
              if(k > 1 && k < (N - 1))
              {
                if ((data_ir_vals[k + jump] > 0) && (data_ir_vals[k+ jump + 1] > 0))
                {
                  //Set Starting Point to before the two NonZero Points
                  truncAddress[i] = (k - 1);
                  break;
                }
                if ((data_ir_vals[k + jump] < 0) && (data_ir_vals[k+ jump +1]) < 0)
                {
                  truncAddress[i] = (k - 1);
                  break;
                }
              }

              //**********************
              //if two nonzero point happen before second value - set First memory location as the memory allocation to start truncating from
              else if(k <= 1 )
              {
                if (data_ir_vals[k + jump] > 0 && data_ir_vals[k + jump +1] > 0)
                {
                  //Set Starting Point to before the two NonZero Points
                  truncAddress[i] = 0;
                  break;
                }
                else if (data_ir_vals[k + jump] < 0 && data_ir_vals[k + jump +1] < 0)
                {
                  //Set Starting Point to before the two NonZero Points
                  truncAddress[i] = 0;
                  break;
                }
              }

              else if(k == (N - 1))
              {
                if (data_ir_vals[k + jump] > 0 && data_ir_vals[k + jump -1] > 0)
                {
                  //Set Starting Point to before the two NonZero Points
                  truncAddress[i] = 0;
                  break;
                }
                else if (data_ir_vals[k + jump] < 0 && data_ir_vals[k + jump -1] < 0)
                {
                  //Set Starting Point to before the two NonZero Points
                  truncAddress[i] = 0;
                  break;
                }

                }
              }
              if(j <(R-1))
              {
                jump += N;
              }
            }
            jump +=N;
          }

      //Test Values by Print
      /*for(i=0; i<M;i++)
      {
        printf("\n\nTruncAddress %i = %i", i, truncAddress[i]);
      }
      */


      printf("\n\n Finding earliest address - common for consistent delay\n\n");

      //**********************
      //Find Minimum Value in TruncAddress
      minimumDelayValue = truncBits;
      for (i=0;i<M;i++)
      {
        if(truncAddress[i] < minimumDelayValue)
        {
          minimumDelayValue = truncAddress[i];
        }
      }

      //printf("\n\nAll working fine before Value Transfer\n\n");
      jump = 0;
      jumpb = 0;
      for(i=0;i<M;i++) //Cycles through recievers
      {
          for (j=0;j<R;j++) //Cycles through Measurements
          {
            for(k=minimumDelayValue;k<(truncBits+minimumDelayValue);k++) //Cycles through truncated Length samples
              {
                data_ir_valsTrunc[(k-minimumDelayValue)+jump] = data_ir_vals[k+jumpb];
              }

              if (j < (R-1))
              {
                jumpb += N;
                jump += truncBits;
              }
          }
      jumpb += N;
      jump += truncBits;
      }

      //Add fade out to last ten samples
      jump = 0;
      for(i=0;i<M;i++) //Cycles through Measurements
      {
          for (j=0;j<R;j++) //Cycles through Recievers
          {
            for(k=0;k<truncBits;k++) //Cycles through Samples
            {
              if(k >= (truncBits - 10) && k < truncBits)
              {
                data_ir_valsTrunc[k+jump] *= ((truncBits - (k+1))/10);
              }

            }
            if (j < (R-1))
            {
              jump += truncBits;
            }
          }
          jump += truncBits;
      }

    //printf("\n\n Fade out added\n\n");
  }//end of if statement when interpol is not needed and number of samples is not equal

      if(interpolNeed == 1 && N == truncBits)
      {
        jump = 0;
        jumpb = 0;
        for(i=0;i<M;i++) //Cycles through Measurements
        {
            for (j=0;j<R;j++) //Cycles through Recievers
            {
              for(k=0;k<truncBits;k++) //Cycles through Samples
              {
                  data_ir_valsTrunc[k+jump] = data_ir_vals[k+jumpb];
                  //printf("\n\ndata_ir_valsTrunc = %f\n\n", data_ir_valsTrunc[k+jump]);
              }

              if (j < (R-1))
              {
                jumpb += N;
                jump += truncBits;
              }
            }
            jumpb += N;
            jump += truncBits;
        }
      }

      if(interpolNeed == 0 && N == truncBits)
      {
        jump = 0;
        for(i=0;i<M;i++) //Cycles through Measurements
        {
            for (j=0;j<R;j++) //Cycles through Receivers
            {
              for(k=0;k<truncBits;k++) //Cycles through Samples
              {
                  data_ir_valsTrunc[k+jump] = data_ir_vals[k+jump];
                  //printf("\n\ndata_ir_valsTrunc = %f\n\n", data_ir_valsTrunc[k+jump]);
              }

              if (j < (R-1))
              {
                jump += truncBits;
              }
            }
            jump += truncBits;
        }
      }

      /*
        jump = 0;
        for(i=0;i<M;i++) //Cycles through recievers
        {
            for (j=0;j<R;j++) //Cycles through Measurements
            {
              for(k=0;k<truncBits;k++) //Cycles through Samples
              {
                  printf("\n\ndata_ir_valsTrunc = %f\n\n", data_ir_valsTrunc[k+jump]);
              }

              if (j < (R-1))
              {
                jump += truncBits;
              }
            }
            jump += truncBits;
        }
        */


/*****************/
/*****************/
//*INTERPOLATION*//
/*****************/
/*****************/
        printf("\n\nchecking distances from measurements to theoretical points for interpolation\n\n");
        //setting up an array to store all the distances from the source positions to our current point
        distance = malloc(((SourcePositionLength[0]*SourcePositionLength[1])/3) * sizeof(double));
        //printf("\n\n Finding closest points\n\n");

        //If interpolation is needed, and the data_ir_values have been truncated
        if(interpolNeed == 1 && data_ir_valsTrunc != NULL)
        {
          jump = 0;
          azimuth = 0;
          zenith = 0;

          //getting the distance between all points in the source, and the MIT data set
          for(i=0; i<710*3; i+=3)
          {
            TestAngleHi = 0;
            TestAngleLo = 0;
            //Equation below equivalent to Radial Distance *sin(Zenith) *cos(azimuth)
            //Simplify by assuming radial distance = 1 for all sources (distance is emulated in the opcode)
            //Calculate the cartesian positions
            xa = cos((CompareSourcePosition[i+1]*M_PI)/180) * cos((CompareSourcePosition[i]*M_PI)/180);
            ya = cos((CompareSourcePosition[i+1]*M_PI)/180) * sin((CompareSourcePosition[i]*M_PI)/180);
            za = sin((CompareSourcePosition[i+1]*M_PI)/180);

                  //Cycle through all source positions
                  for(k=0;k<(SourcePositionLength[0]*SourcePositionLength[1]);k+=3)
                  {
                    azimuth = (SourcePosition[k]*M_PI)/180;
                    zenith = (SourcePosition[k+1]*M_PI)/180;
                    //fprintf(ValDistTest,"\n\nSOURCE POSITION = %f %f %f\n\n",SourcePosition[k],SourcePosition[k+1],SourcePosition[k+2]);
                    //fprintf(ValDistTest,"\n\nRotation angle= %f polar angle = %f\n\n",azimuth,zenith);

                    //Simplify by assuming radial distance = 1
                    x = cos(zenith) * cos(azimuth);
                    y = cos(zenith) * sin(azimuth);
                    z = sin(zenith);

                    //this should fill a M size array with the distance from each source point to the MIT Data point
                    //with the location of the distance in the array equal to its measurement "number"

                    distance[(k/3)] = sqrt(SQUARE(xa-x) + SQUARE(ya-y) + SQUARE(za-z));
                    //fprintf(ValDistTest,"\n\nMIT A = %f MIT Z = %f MIT RD = %f, \nsource A = %f source Z = %f source RD = %f",CompareSourcePosition[i],CompareSourcePosition[i+1],CompareSourcePosition[i+2],SourcePosition[k],SourcePosition[k+1],SourcePosition[k+2]);
                    //fprintf(ValDistTest,"\n\nXA = %f YA = %f ZA = %f, x = %f y = %f z = %f",xa,ya,za,x,y,z);
                    //fprintf(ValDistTest,"\n\nDistance from MIT %i for measurement %i = %f\n\n",i/3,k/3,distance[k/3]);
                  }
            //Set the "min array" to be the first four values of the array "distance"
            min[0] = distance[0];
            min[1] = distance[1];
            min[2] = distance[2];
            min[3] = distance[3];

            location[0] = 0;
            location[1] = 1;
            location[2] = 2;
            location[3] = 3;

            //Sort these into order of minimum value - with a the location array storing the orignal location of each point in the distance array
            for (l = 0; l < 4; l++)
            {
              for (m = l + 1; m < 4; m++)
              {
                if (min[l] > min[m])
                {
                  temp =  min[l];
                  min[l] = min[m];
                  //location[l] = m;
                  min[m] = temp;
                  //location[m] = l;
                }
              }
            }
            //THiS COULD BE AN ISSUE WITH TWO POINTS WITH Exact SAME DISTANCE in the first four measurements


            for(l=0; l <4; l++)
            {
              if (min[l]==distance[0])
              {
                location[l] = 0;
              }
              else if(min[l]==distance[1])
              {
                location[l] = 1;
              }
              else if(min[l]==distance[2])
              {
                location[l] = 2;
              }
              else if(min[l]==distance[3])
              {
                location[l] = 3;
              }
            }
            //Test By Print
            /*printf("\n\nPoints closest to MIT measurement %i are stored at\n",(i/3));
            for (p=0; p<4; p++)
            {
                 printf("\nFirst four distances - Location %i = %i\n",p,location[p]);
                 printf("\nwith a distance of %f\n\n",min[p]);
            }
            */


            //Compare the rest of the values to these four storing the value, and its original location into arrays
            for (l = 4; l <(SourcePositionLength[0]*SourcePositionLength[1])/3; l++)
            {
              TestAngleHi = CompareSourcePosition[(i)] + 45;
              if (TestAngleHi > 360)
              {
                TestAngleHi -= 360;
              }
              TestAngleLo = CompareSourcePosition[i] - 45;
              if (TestAngleLo < 0)
              {
                TestAngleLo += 360;
              }
              SPAngle = SourcePosition[(l*3)];

              if(CompareSourcePosition[i] < 45 || CompareSourcePosition[i] > 315)
              {
                if((SPAngle >= 0 && SPAngle <= TestAngleHi) || (SPAngle >= TestAngleLo && SPAngle < 360))
                {
                  ignore = 0;
                }
                else
                {
                  ignore = 1;
                }
              }
              else
              {
                if((SPAngle <= TestAngleHi) && (SPAngle >= TestAngleLo))
                {
                  ignore = 0;
                }
                else
                {
                  ignore = 1;
                }
              }

              if (distance[l] < min[0] && ignore == 0)
              {
                  min[3] = min[2];
                  location[3] = location[2];
                  min[2] = min[1];
                  location[2] = location[1];
                  min[1] = min[0];
                  location[1] = location[0];
                  min[0] = distance[l];
                  location[0] = l;
              }
              else if (distance[l] < min[1] && distance[l] >= min[0] && ignore == 0)
              {
                  min[3] = min[2];
                  location[3] = location[2];
                  min[2] = min[1];
                  location[2] = location[1];
                  min[1] = distance[l];
                  location[1] = l;
              }
              else if (distance[l] < min[2] && distance[l] >= min[1] && ignore == 0)
              {
                  min[3] = min[2];
                  location[3] = location[2];
                  min[2] = distance[l];
                  location[2] = l;
              }
              else if (distance[l] < min[3] && distance[l] >= min[2] && ignore == 0)
              {
                   min[3] = distance[l];
                   location[3] = l;
              }
             }

               //Test By Print
               /*printf("\n\nPoints closest to MIT measurement %i are stored at\n",(i/3));
               for (p=0; p<4; p++)
               {
                    printf("\nLocation %i = %i\n",p,location[p]);
                    printf("\nwith a distance of %f\n\n",min[p]);

                }
                */

                //Calculate the total distance from the current point
                totaldistance = min[0] + min [1] + min [2] + min [3];

                weightA = (totaldistance/min[0]);
                weightB = (totaldistance/min[1]);
                weightC = (totaldistance/min[2]);
                weightD = (totaldistance/min[3]);

                totalweight = weightA + weightB + weightC + weightD;

                //This is used to give the weight based on distance
                weightA = roundf((weightA/totalweight)* 100)/100;
                weightB = roundf((weightB/totalweight)* 100)/100;
                weightC = roundf((weightC/totalweight)* 100)/100;
                weightD = roundf((weightD/totalweight)* 100)/100;

                if(min[0] == 0)
                {
                  weightA = 1;
                  weightB = 0;
                  weightC = 0;
                  weightD = 0;
                }

                if ((weightA + weightB + weightC + weightD) != 1)
                {
                  weightA = 1 - (weightB + weightC + weightD);
                }

                //Test by Print
                //printf("\n\nMeasurement %i= A %f B %f C %f D %f\n\n",i/3, weightA,weightB, weightC, weightD);

                //Get the values of each closest sample (with at least 2 recievers)
                //gives values for first two recievers
                if(R >1)
                {
                  for(l=0;l<truncBits;l++)
                  {
                    interpolAl[l] = data_ir_valsTrunc[(location[0]*(R*truncBits))+truncBits+l];
                    interpolAr[l] = data_ir_valsTrunc[(location[0]*(R*truncBits)) + l];

                    interpolBl[l] = data_ir_valsTrunc[(location[1]*(R*truncBits))+truncBits +l];
                    interpolBr[l] = data_ir_valsTrunc[(location[1]*(R*truncBits))+ l] ;

                    interpolCl[l] = data_ir_valsTrunc[(location[2]*(R*truncBits))+truncBits+l];
                    interpolCr[l] = data_ir_valsTrunc[(location[2]*(R*truncBits)) + l];

                    interpolDl[l] = data_ir_valsTrunc[(location[3]*(R*truncBits))+truncBits+l];
                    interpolDr[l] = data_ir_valsTrunc[(location[3]*(R*truncBits)) + l];
                  }
                }

                //Get the values of each closest sample (with only 1recievers)
                //copies value for first reciever, to the second reciever in the output.
                if(R < 2)
                {
                  printf("\n\nWarning ONLY ONE RECIEVER - BINAURAL WILL BE COMPROMISED\n");
                  for(l=0;l<truncBits;l++)
                  {
                    interpolAl[l] = data_ir_valsTrunc[(location[0]*truncBits)+l];
                    interpolAr[l] = data_ir_valsTrunc[(location[0]*truncBits)+l];

                    interpolBl[l] = data_ir_valsTrunc[(location[1]*truncBits)+l];
                    interpolBr[l] = data_ir_valsTrunc[(location[1]*truncBits)+l];

                    interpolCl[l] = data_ir_valsTrunc[(location[2]*truncBits)+l];
                    interpolCr[l] = data_ir_valsTrunc[(location[2]*truncBits)+l];

                    interpolDl[l] = data_ir_valsTrunc[(location[3]*truncBits)+l];
                    interpolDr[l] = data_ir_valsTrunc[(location[3]*truncBits)+l];

                  }
                }

                //each of the closest values gets translated into Freq spectra
                //prior to interpolation
                fftw_execute(SampAplanl);
                fftw_execute(SampAplanr);
                fftw_execute(SampBplanl);
                fftw_execute(SampBplanr);
                fftw_execute(SampCplanl);
                fftw_execute(SampCplanr);
                fftw_execute(SampDplanl);
                fftw_execute(SampDplanr);

                //Test by Print
                /*
                for (l=0; l<truncBits; l++)
                {
                  printf("\n\nSample a left fft results %f\n\n", fftSampAl[l]);
                  printf("\n\nSample a right fft results %f\n\n", fftSampAr[l]);
                  printf("\n\nSample b left fft results %f\n\n", fftSampBl[l]);
                  printf("\n\nSample b right fft results %f\n\n", fftSampBr[l]);
                  printf("\n\nSample c left fft results %f\n\n", fftSampCl[l]);
                  printf("\n\nSample c right fft results %f\n\n", fftSampCr[l]);
                  printf("\n\nSample d left fft results %f\n\n", fftSampDl[l]);
                  printf("\n\nSample d right fft results %f\n\n", fftSampDr[l]);
                }
                */

                /* Set first two values in array to value at 0Hz and nyq */
                polarAl[0] = fftSampAl[0];
                polarAl[1] = fftSampAl[truncBits/2];
                polarAr[0] = fftSampAr[0];
                polarAr[1] = fftSampAr[truncBits/2];

                polarBl[0] = fftSampBl[0];
                polarBl[1] = fftSampBl[truncBits/2];
                polarBr[0] = fftSampBr[0];
                polarBr[1] = fftSampBr[truncBits/2];

                polarCl[0] = fftSampCl[0];
                polarCl[1] = fftSampCl[truncBits/2];
                polarCr[0] = fftSampCr[0];
                polarCr[1] = fftSampCr[truncBits/2];

                polarDl[0] = fftSampDl[0];
                polarDl[1] = fftSampDl[truncBits/2];
                polarDr[0] = fftSampDr[0];
                polarDr[1] = fftSampDr[truncBits/2];

                /* fill the array with calculated mag/phase: polar*/
                for(l = 2, k = 1; l < truncBits; k++, l += 2)
                {
                  polarAl[l] = sqrt(SQUARE(fftSampAl[k]) + SQUARE(fftSampAl[truncBits - k]));
                  polarAl[l+1] = atan2(fftSampAl[truncBits-k],fftSampAl[k]);
                  polarAr[l] = sqrt(SQUARE(fftSampAr[k]) + SQUARE(fftSampAr[truncBits - k]));
                  polarAr[l+1] = atan2(fftSampAr[truncBits-k],fftSampAr[k]);

                  polarBl[l] = sqrt(SQUARE(fftSampBl[k]) + SQUARE(fftSampBl[truncBits - k]));
                  polarBl[l+1] = atan2(fftSampBl[truncBits-k],fftSampBl[k]);
                  polarBr[l] = sqrt(SQUARE(fftSampBr[k]) + SQUARE(fftSampBr[truncBits - k]));
                  polarBr[l+1] = atan2(fftSampBr[truncBits-k],fftSampBr[k]);

                  polarCl[l] = sqrt(SQUARE(fftSampCl[k]) + SQUARE(fftSampCl[truncBits - k]));
                  polarCl[l+1] = atan2(fftSampCl[truncBits-k],fftSampCl[k]);
                  polarCr[l] = sqrt(SQUARE(fftSampCr[k]) + SQUARE(fftSampCr[truncBits - k]));
                  polarCr[l+1] = atan2(fftSampCr[truncBits-k],fftSampCr[k]);

                  polarDl[l] = sqrt(SQUARE(fftSampDl[k]) + SQUARE(fftSampDl[truncBits - k]));
                  polarDl[l+1] = atan2(fftSampDl[truncBits-k],fftSampDl[k]);
                  polarDr[l] = sqrt(SQUARE(fftSampDr[k]) + SQUARE(fftSampDr[truncBits - k]));
                  polarDr[l+1] = atan2(fftSampDr[truncBits-k],fftSampDr[k]);

                  /*printf("\n\npolarAl = %f\n\n", polarAl[l]);
                  printf("\n\npolarAr = %f\n\n", polarAr[l]);
                  printf("\n\npolarBl = %f\n\n", polarBl[l]);
                  printf("\n\npolarBr = %f\n\n", polarBr[l]);
                  printf("\n\npolarCl = %f\n\n", polarCl[l]);
                  printf("\n\npolarCr = %f\n\n", polarCr[l]);
                  printf("\n\npolarDl = %f\n\n", polarDl[l]);
                  printf("\n\npolarDr = %f\n\n", polarDr[l]);
                  */
                }

                //0Hz and Nyquist were never transformed into polar co-ordinates - pos/negetive value dictates phase
                //Use absolute values of magnitue (to negate the phase) then check nearest point to add this back in
                HRTFfreql[0] = fabs(weightA*polarAl[0]) + fabs(weightB*polarBl[0]) + fabs(weightC*polarCl[0]) + fabs(weightD*polarDl[0]);
                HRTFfreqr[0] = fabs(weightA*polarAr[0]) + fabs(weightB*polarBr[0]) + fabs(weightC*polarCr[0]) + fabs(weightD*polarDr[0]);

                if (polarAl[0] < 0)
                {
                  HRTFfreql[0] = -HRTFfreql[0];
                }
                else
                {
                  HRTFfreql[0] = fabs(HRTFfreql[0]);
                }
                if (polarAr[0] < 0)
                {
                  HRTFfreqr[0] = -HRTFfreqr[0];
                }
                else
                {
                  HRTFfreqr[0] = fabs(HRTFfreqr[0]);
                }

                HRTFfreql[1] = fabs(weightA*polarAl[1]) + fabs(weightB*polarBl[1]) + fabs(weightC*polarCl[1]) + fabs(weightD*polarDl[1]);
                HRTFfreqr[1] = fabs(weightA*polarAr[1]) + fabs(weightB*polarBr[1]) + fabs(weightC*polarCr[1]) + fabs(weightD*polarDr[1]);

                if (polarAl[1] < 0)
                {
                  HRTFfreql[1] = -HRTFfreql[1];
                }
                else
                {
                  HRTFfreql[1] = fabs(HRTFfreql[1]);
                }
                if (polarAr[1] < 0)
                {
                  HRTFfreqr[1] = -HRTFfreqr[1];
                }
                else
                {
                  HRTFfreqr[1] = fabs(HRTFfreqr[1]);
                }

                //Interpolation - weighs of magnitude of each sample
                for(l=2; l<truncBits; l+=2)
                {
                  HRTFfreql[l] = (weightA*polarAl[l]) + (weightB*polarBl[l]) + (weightC*polarCl[l]) + (weightD*polarDl[l]);
                  HRTFfreqr[l] = (weightA*polarAr[l]) + (weightB*polarBr[l]) + (weightC*polarCr[l]) + (weightD*polarDr[l]);

                  //printf("\n\nHRTFfreql Mag %i= %f\n\n",i/3, HRTFfreql[l]);
                  //printf("\n\nHRTFfreqr Mag %i= %f\n\n",i/3, HRTFfreqr[l]);
                }
                //Store the phase values of the closest element
                for(l=3; l<truncBits; l+=2)
                {
                  HRTFfreql[l] = polarAl[l];
                  HRTFfreqr[l] = polarAr[l];
                  //printf("\n\nHRTFfreql Phase %i = %f\n\n",l, HRTFfreql[l]);
                  //printf("\n\nHRTFfreqr PHase %i= %f\n\n",l, HRTFfreqr[l]);
                }

                //HRTFfreql & HRTFfreqr are the value of the dataset for storing
                for(l=0;l<truncBits; l++)
                {
                data_HRTF_Outputl[jump+l] = (float) HRTFfreql[l];
                data_HRTF_Outputr[jump+l] = (float) HRTFfreqr[l];
                //printf("\n\ndata_HRTF_Output LEFT = %f\n\n", data_HRTF_Outputl[jump+l]);
                //printf("\n\ndata_HRTF_Output RIGHT = %f\n\n", data_HRTF_Outputr[jump + l]);
                }
                jump += truncBits;
          }

          /*ValOutTest = fopen("Testdata.txt", "wb");
          for(i=0;i<710*128;i++)
          {
            fprintf(ValOutTest,"\n\ndata_HRTF_Output left %i = %f\n\n",i,data_HRTF_Outputl[i]);
            fprintf(ValOutTest,"\n\ndata_HRTF_Output right %i = %f\n\n",i,data_HRTF_Outputr[i]);
          }
          */
        }


        //If interpolation is not needed, and the data_ir_values have been truncated
        // HRIR impulses still need to be fft'd to Frequency Spectra
        if(interpolNeed == 0 && data_ir_valsTrunc != NULL)
        {
          printf("interpol not needed, FFTing data");
          jump = 0;
          jumpb =0;
          //Need to change all values to frequency spectra
          for(i=0; i<710; i++)
          {
                //Get the values of each sample (with at least 2 recievers)
                //gives values for first two recievers
                if(R >1)
                {
                  for(l=0;l<truncBits;l++)
                  {
                    interpolAl[l] = data_ir_valsTrunc[jumpb + truncBits +l];
                    interpolAr[l] = data_ir_valsTrunc[jumpb  + l];

                    //Test by Print - Fault find
                    //printf("\n\nSample a left %f\n\n", interpolAl[l]);
                    //printf("\n\nSample a right %f\n\n", interpolAr[l]);
                  }
                  jumpb +=truncBits;
                }

                //Get the values of each sample (with only 1recievers)
                //copies value for first reciever, to the second reciever in the output.
                if(R < 2)
                {
                  for(l=0;l<truncBits;l++)
                  {
                    interpolAl[l] = data_ir_valsTrunc[l + jumpb];
                    interpolAr[l] = data_ir_valsTrunc[l + jumpb];
                  }
                }

                //each of the closest values gets translated into Freq spectra
                //prior to interpolation
                fftw_execute(SampAplanl);
                fftw_execute(SampAplanr);


                /* Set first two values in array to value at 0Hz and nyq */
                interpolAl[0] = fftSampAl[0];
                interpolAl[1] = fftSampAl[truncBits/2];
                interpolAr[0] = fftSampAr[0];
                interpolAr[1] = fftSampAr[truncBits/2];


                /* fill the array with calculated mag/phase: polar*/
                for(l = 2, k = 1; l < truncBits; k++, l += 2)
                {
                  interpolAl[l] = sqrt(SQUARE(fftSampAl[k]) + SQUARE(fftSampAl[truncBits - k]));
                  interpolAl[l+1] = atan2(fftSampAl[truncBits-k],fftSampAl[k]);
                  interpolAr[l] = sqrt(SQUARE(fftSampAr[k]) + SQUARE(fftSampAr[truncBits - k]));
                  interpolAr[l+1] = atan2(fftSampAr[truncBits-k],fftSampAr[k]);
                }

                //interpolAl & interpolAr are the value of the dataset for storing
                for(l=0;l<truncBits; l++)
                {
                data_HRTF_Outputl[jump+l] = (float) interpolAl[l];
                data_HRTF_Outputr[jump+l] = (float) interpolAr[l];
                }

                jumpb += truncBits;
                jump += truncBits;
          }
        }

        /*
        ValOutTest = fopen("Testdata.txt", "wb");
        for(i=0;i<710*128;i++)
        {
          fprintf(ValOutTest,"\n\ndata_HRTF_Output l%i = %f\n\n",i,data_HRTF_Outputl[i]);
          fprintf(ValOutTest,"\n\ndata_HRTF_Output r%i = %f\n\n",i,data_HRTF_Outputr[i]);
        }
        */

        //Write data_HRTF_Output to a new file
        stringwrite(data_SampleRate);
        foutl = fopen(filenamel, "wb");
      	foutr = fopen(filenamer, "wb");

        jump = 0;
        noofsamplesoutput = 0;
        pdata_HRTF_Outputl = &data_HRTF_Outputl[0];
        pdata_HRTF_Outputr = &data_HRTF_Outputr[0];
        for(i=0;i<710*3; i+=3)
        {
          if (CompareSourcePosition[i] <= 180)
          {
            noofsamplesoutput++;
            //printf("\n\nSourcePosition = %f\n\n",CompareSourcePosition[i]);
            fwrite(pdata_HRTF_Outputl, sizeof(float), truncBits, foutl);
    		    fwrite(pdata_HRTF_Outputr, sizeof(float), truncBits, foutr);
          }
          jump += truncBits;
          pdata_HRTF_Outputl = &data_HRTF_Outputl[jump];
          pdata_HRTF_Outputr = &data_HRTF_Outputr[jump];
        }

        printf("\n\nNo of samples printed to the .dat files %i\n\n", noofsamplesoutput);
        //samplenumber = 0;
        //for(i=0;i<noofsamplesoutput*128;i++)
        //{
        //  printf("\n\ndata_HRTF_Output l%i = %f\n\n",samplenumber,data_HRTF_Outputl[i]);
        //  printf("\n\ndata_HRTF_Output r%i = %f\n\n",samplenumber,data_HRTF_Outputr[i]);
        //}


      free(data_ir_vals);
      free(ListenerPosition);
      free(SourcePosition);
      free(SpherePosition);
      free(EmitterPosition);
      free(data_ir_valsTrunc);
      free(truncAddress);
      free(distance);
      fftw_destroy_plan(SampAplanl);
      fftw_destroy_plan(SampAplanr);
      fftw_destroy_plan(SampBplanl);
      fftw_destroy_plan(SampBplanr);
      fftw_destroy_plan(SampCplanl);
      fftw_destroy_plan(SampCplanr);
      fftw_destroy_plan(SampDplanl);
      fftw_destroy_plan(SampDplanr);
      fclose(foutl);
    	fclose(foutr);
  }
  else
  {
    printf("\n\nFile Does not exist in the current folder, program cannot be run\n\n");
  }


  printf("\n\nFINISHED!!!!\n\n");

}
