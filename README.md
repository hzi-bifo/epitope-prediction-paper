#  “EpitopeVec:   LinearEpitope   Prediction   Using   DeepProtein   Sequence   Embeddings“
Data,sscripts and results for the EpitopeVec article by Bahai et al., in review.

  The epitope-prediction software is available at https://github.com/hzi-bifo/epitope-prediction

## Requirements

* **```Python 3```** with the following packages:
    * **numpy 1.17.1**
    
 ## Usage
 
 * Clone this repository:
  ```
  git clone https://github.com/hzi-bifo/epitope-prediction-paper
  ```

Run the read_prediction.py file with the pickle file as an argument. Each folder contains the results for each dataset with their respective pickle file.
       
       python3 read_prediction.py ./ABCPred/abcpred.pickle
      
The results with accuracy, precision, recall, ROC metrics etc. will be displayed.
