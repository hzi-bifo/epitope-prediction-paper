#  EpitopeVec:   LinearEpitope   Prediction   Using   DeepProtein   Sequence   Embeddings
Data, scripts and results for the EpitopeVec article by Bahai et al., Bioinformatics, 2021.

  The epitope-prediction software is available at https://github.com/hzi-bifo/epitope-prediction

## Requirements

* **```Python 3```** with the following packages:
    * **numpy 1.17.1**
    * **scipy 1.4.1**
    * **matplotlib 3.1.3**
    * **sklearn 0.22.1**
    * **pydpi 1.0**
    * **biopython 1.71.0**
    * **tqdm 4.15.0**
    * **gensim 3.8.3**
    
   
  If these are not installed, you can install them with ``` pip ```. 
    ```
   pip3 install -r ./requirement/requirements.txt
   ```
   
  Additionally, **pydpi 1.0** from ```pip``` might be incompatible with **Python 3**. Please install the **pydpi** package from the provided ```pydpi.tar.gz``` file.
    ```
    pip3 install pydpi.tar.gz
    ```
   
 * Binary file for ProtVec representation of proteins can be downloaded using the following command in the ```protvec``` directory:
 
 ```
 cd protvec
 wget http://deepbio.info/embedding_repo/sp_sequences_4mers_vec.txt
 wget http://deepbio.info/embedding_repo/sp_sequences_4mers_vec.txt.bin -O sp_sequences_4mers_vec.bin
 ```
    
## Usage
 
* Clone this repository:
  ```
  git clone https://github.com/hzi-bifo/epitope-prediction-paper
  ```

* To train a new machine learning model, run the training file with name of the dataset you want to train on. The datasets are inside the retraining folder (bcpreds, ibce-el, lbtope and viral). eg:
  ```
   python3 retrain.py bcpreds
  ```
  
  Use ```bcpreds``` for training on the BCPreds dataset.
  Use ```lbtope``` for training on the LBTope dataset.
  Use ```ibce-el``` for training on the iBCE-EL training dataset.


## Input

* For training a new model, two files containing a list of confirmed positive and negative epitopes are needed. These can be .txt files with each line containing a peptide. eg: In the **./retraining/bcpreds/** folder **pos.txt** contains a list of petides which are epitopes and **neg.txt** contains non-epitopes.  
For training domain-specific models, all the epitope petides should also be from the specific domain. eg: If one wants a viral-specific model, only include epitopes derived from viral proteins.


## Output

* Training a new model will create a pickle file in the **/retraining/input dataset** folder. The **modelname.pickle** is the newly trained model which can be used with the EpitopeVec software (https://github.com/hzi-bifo/epitope-prediction) for testing.

## Testing

* To test the performance of the trained models, please use the **test.py** file.

  ```
  python3 test.py model_pickle_file peptidefile

  ```

Here,

**model_pickle_file** is the _.pickle_ file you trained previosly.
**peptidefile** is a file with two columns. The first column is the peptide sequences and the second is the target (1 for epitope and 0 for non-epitope). See **testing** directory for a list of example files. The two columns should be tab-separated.
Please provide the full name and location of the model\_pickle file and the full name and location of the peptidefile. For example, if you want to test the model trained on the iBCE-EL datatset (svm-ibce-el.pcikle) on the ABCPred16 dataset run 
```
python3 test.py ./retraining/ibce-el/svm-ibce-el.pickle ./testing/abcpred16.txt 

```

## Results

The becnhmarking results of various methods are inside the **results** folder. 



