import os,sys

## ROUTINES FOR AUTOMATIC PATH SETTING
def clean_path(cpath):
  npath = cpath.replace('/./','/').replace('//','/')
  while (cpath != npath):
    cpath = npath
    npath = cpath.replace('/./','/').replace('//','/')
  return npath
if (len(sys.argv[0]) > 0)  and os.path.isdir(os.getcwd()+'/'+os.path.dirname(sys.argv[0])):
  def_ibis_path = os.path.dirname(os.getcwd()+'/'+sys.argv[0])+'/'
elif (len(sys.argv[0]) > 0)  and os.path.isdir(os.path.dirname(sys.argv[0])):
  def_ibis_path = os.path.dirname(sys.argv[0])+'/'
else:
  def_ibis_path = os.getcwd()+'/'
def_ibis_path=clean_path(def_ibis_path)

## ALTERNATIVE MANUAL PATH SETTING
#def_ibis_path = '/mnt/solexa/bin/'

def_soap_path=def_ibis_path+'soap_1.11_patched/soap'
def_soap_ref='/mnt/solexa/Genomes/phiX/phiX174.fa'

def_bowtie_path='PATH_TO_BOWTIE/bowtie'
def_bowtie_ref='/mnt/solexa/Genomes/phiX/bowtie'

def_extract_dataset = def_ibis_path+"createTrainingSeqs.py" ##SOAP V1
def_align_path = def_soap_path
def_align_ref = def_soap_ref
#def_extract_dataset = def_ibis_path+"createTrainingSeqs_bowtie.py" ## BOWTIE 0.12.4 (tested)
#def_align_path = def_bowtie_path
#def_align_ref = def_bowtie_ref


#TRAINING PARAMS
maxTrainTime=1000;
def_training = def_ibis_path+"trainingSeqs2SVMlightPerCycle.py";
plot_recal   = def_ibis_path+"plot_recal.R";

def_svmlight_path =def_ibis_path+'libocas_v093/';
def_predictor_path =def_ibis_path+'predictor/';

def_liblinear_path=def_ibis_path+'liblinear-1.8mod/';

def_svm_train = def_svmlight_path+"msvmocas -v 0 -t "+str(maxTrainTime);
def_log_train = def_liblinear_path+"train -q -B 0 -e 0.001 -s 2  ";

def_svm_test = def_svmlight_path+'linclass';

#QUALITY SCORE recalibration
def_MAXQUALSCORE=40;
def_MINQUALSCORE=1;
def_recalib_path        = def_ibis_path+'qualRecal/'
def_estimateError_bin   = def_recalib_path+'estimaErrorLookupMax';
def_estimateError       = def_estimateError_bin+" 50 75 45 "; #[minimum number of datapoints] [maximum number of datapoints] [return a datapoint after having reached that value]
def_linearRegression_bin= def_recalib_path+'findBestk';
def_linearRegression    = def_linearRegression_bin+" ";
maxTrainTimeBinary=1000;

def_svm_trainbinary_bin = def_liblinear_path+"train";
def_svm_trainbinary = def_svm_trainbinary_bin+" -q -s 0 ";
def_svm_testbinary_bin  = def_svmlight_path+'linclass';
def_svm_testbinary      = def_svm_testbinary_bin+' -t 1';


#MASKING PARAMS
def_soapmask = def_ibis_path+"soap2masked.py"
def_percentageForMask=10;
def_qualityForMask   =0;



def_svm_prediction = def_predictor_path+'svm_libocas_classify';
def_prediction_py = def_ibis_path+"runPredictionOnly.py"
def_temp='/tmp/'




