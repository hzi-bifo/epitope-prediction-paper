#  EpitopeVec:   LinearEpitope   Prediction   Using   DeepProtein   Sequence   Embeddings
Data,scripts and results for the EpitopeVec article by Bahai et al., in review.

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
    * **gensim 2.3.0**
    
   
  If these are not installed, you can install them with ``` pip ```. 
    ```
   pip3 install -r ./requirement/requirements.txt
   ```
   
  Additionally, **pydpi 1.0** from ```pip``` might be incompatible with **Python 3**. Please install the **pydpi** package from the provided ```pydpi.tar.gz``` file.
    ```
    pip3 install pydpi.tar.gz
    ```
   
 * Binary file for ProtVec representation of proteins. It can be downloaded from https://drive.google.com/drive/folders/1Y3PgFzUwCKu7Jzn0XSwd7UaBmDkXpTBZ?usp=sharing and then it should be put in the ```protvec``` folder.
    
## Usage
 
* Clone this repository:
  ```
  git clone https://github.com/hzi-bifo/epitope-prediction-paper
  ```

* To train a new machine learning model, run the training file with the positive(epitope) and negative (non-epitope set)
  ```
   python3 train.py -p positiveset.txt -n negativeset.txt
  ```

* To see the results for different datasets, run the **read_prediction** file with the corresponding **.pip** file as input argument. eg:
  ```
  python3 read_prediction.py ./ABCPred/abcpred.pickle
  ```

## Input

* For training a new model, two files containing a list of confirmed positive and negative epitopes are needed. These can be .txt files with each line containing a peptide. eg: In the **ABCPred** folder **abcpred_pos_20mer.txt** contains a list of petides which are epitopes and **abcpred_neg_20mer.txt** contains non-epitopes.  
For training domain-specific models, all the epitope petides should also be from the specific domain. eg: If one wants a viral-specific model, only include epitopes derived from viral proteins.

* Each dataset contains a **.pickle** file. These pickle files will reproduce the results mentioned in the paper on running them with the **read_prediction** file.

## Output

* Training a new model will create two pickle files in the **model** folder. The **model.pickle** is the newly trained model which can be used with the EpitopeVec software (https://github.com/hzi-bifo/epitope-prediction) for testing and **scaling.pickle** for scaling the new data being tested.

* Running **read_prediction** file will dipslay the name of the metric and the corresponding values for the dataset corresponding to the input file.
