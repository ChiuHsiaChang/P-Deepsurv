import numpy as np
import pandas as pd
import time
from sklearn.model_selection import train_test_split
data = pd.read_csv('/home/kyro_zhang/ZQX/data/final_data/clinical_gene1000.csv')
data=data.drop(columns=['id'])
train_data,test_data=train_test_split(data,test_size=0.2, random_state=4)
for i in [test_data]:
    i.index = range(i.shape[0])


for i in [train_data]:
    i.index = range(i.shape[0])

from tfdeepsurv.datasets import survival_stats
colname_e = 'fstatus'
colname_t = 'futime'
survival_stats(train_data, t_col=colname_t, e_col=colname_e, plot=False)
survival_stats(test_data, t_col=colname_t, e_col=colname_e, plot=False)
from tfdeepsurv.datasets import survival_df
surv_train = survival_df(train_data, t_col=colname_t, e_col=colname_e, label_col="Y")
surv_test = survival_df(test_data, t_col=colname_t, e_col=colname_e, label_col="Y")
surv_train.head()
from tfdeepsurv import dsnn

#=========================================


def nn(e,f,g,h,_train):
    hidden_layers_nodes = [56,28,14,6,3,1]
    nn_config = {
        #-------0.08 0.03 0.001
        "learning_rate": e,
        "learning_rate_decay": f,
        "activation": 'relu', 
        "L1_reg": g, 
        "L2_reg": h, 
        "optimizer": 'adam',
        "dropout_keep_prob": 1.0,
        "seed": 1
    }
    model = dsnn(
        112, 
        hidden_layers_nodes,
        nn_config
    )
    model.build_graph()
    Y_col = ["Y"]
    X_cols = [c for c in _train.columns if c not in Y_col]
    watch_list = model.train(
        _train[X_cols], _train[Y_col],
        num_steps=900,
        num_skip_steps=30,
        plot=False
    )
    evals_ = model.evals(_train[X_cols], _train[Y_col])
    return evals_


#main
input_nodes = len(surv_train.columns) - 1
hidden_layers_nodes = [5,20,110,1]
nn_config = {
    #-------0.08 0.03 0.001
    "learning_rate": 0.001,
    "learning_rate_decay": 1.0,
    "activation": 'relu', 
    "L1_reg": 5e-06, 
    "L2_reg": 0.05, 
    "optimizer": 'adam',
    "dropout_keep_prob": 1.0,
    "seed": 1
}

model = dsnn(
    input_nodes, 
    hidden_layers_nodes,
    nn_config
)

model.build_graph()

Y_col = ["Y"]
X_cols = [c for c in surv_train.columns if c not in Y_col]
T1 = time.perf_counter()
watch_list = model.train(
    surv_train[X_cols], surv_train[Y_col],
    num_steps=900,
    num_skip_steps=30,
    plot=True
)
T2 =time.perf_counter()
print("CI:", model.evals(surv_train[X_cols], surv_train[Y_col]))
print('cost time:%sms' % ((T2 - T1)*1000))
# model.summary()
#  /anaconda3/envs/R_env/TFDeepSurv/tfdeepsurv
model.close_session()

exit()