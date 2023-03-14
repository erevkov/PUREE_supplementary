import numpy as np
from sklearn.linear_model import LinearRegression

# logit regression definition (credit: https://stackoverflow.com/questions/44234682/how-to-use-sklearn-when-target-variable-is-a-proportion, MB-F)
class LogitRegression(LinearRegression):

    def fit(self, x, p):
        p = np.asarray(p)
        y = np.log(p / (1 - p))
        return super().fit(x, y)

    def predict(self, x):
        y = super().predict(x)
        return 1 / (np.exp(-y) + 1)
    
import tensorflow as tf
from tensorflow import keras

from tensorflow.keras import regularizers
from tensorflow.keras.layers import Dense, Dropout
from tensorflow.keras.models import Sequential
from tensorflow.keras.optimizers import Adam

# fcNN definition
def fcNN(num_features=60000, learning_rate=0.05):
        
        model = Sequential()
        
        if (num_features >= 3000):
            scaling_factor = 10 # decrease substantially if the amount of features is too high
            rounding_factor = 3 # factor to round down
        else:
            scaling_factor = 2
            rounding_factor = 1

        units_init = np.round(num_features, -rounding_factor) // scaling_factor

        # input layer 
        model.add(
            Dense(
                units=units_init,
                activation='relu',
                kernel_regularizer=regularizers.l2(0.01),
                input_dim=num_features))

        model.add(Dropout(rate=0.1))
        
        additional_layers = (num_features // 30000) + 1
        
        # adding intermediate hidden layers
        for i in range(additional_layers):
            units_scaling = scaling_factor ** (2*(i+1))
            model.add(
                Dense(
                    units=(units_init // units_scaling),
                    activation='relu',
                )
            )

        # output layer
        model.add(Dense(1, activation='sigmoid'))

        model.compile(
            optimizer=Adam(learning_rate/1000),
            loss='mean_squared_error')
        
        return model
    
    
### models definitions ###

from sklearn.experimental import enable_halving_search_cv
from sklearn.model_selection import HalvingGridSearchCV
from sklearn.linear_model import LinearRegression, ElasticNetCV, LassoCV, RidgeCV,
from sklearn.svm import NuSVR
from sklearn.ensemble import GradientBoostingRegressor

# Linear Regression #
linreg = LinearRegression()

# Logit Reg #
logitreg = LogitRegression()

# Lasso #
lasso = LassoCV(n_jobs=n_cpu-1, max_iter=100000, cv=5, 
                tol=0.001, selection='random', verbose=True)
# Elastic Net #
elnet = ElasticNetCV(n_jobs=n_cpu-1, max_iter=100000, cv=5,
                     tol=0.001, selection='random', verbose=True)
# Ridge #
ridge = RidgeCV()

# nuSVR #
nusvr = HalvingGridSearchCV(estimator=NuSVR(gamma='auto', kernel='rbf', cache_size=1000), 
                param_grid={'nu':[0.001, 0.01, 0.1, 0.5, 1], 
                            'C':[0.001, 0.01, 0.1, 1, 10]},
                 cv=5,
                 factor=3,
                 scoring='neg_mean_squared_error',
                 n_jobs=n_cpu-1,
                 verbose=3,
                 error_score='raise'
                )

# Gradient Boosting Regressor 
gradboost = HalvingGridSearchCV(estimator=GradientBoostingRegressor(), 
                    param_grid={'max_depth':np.arange(1, 101, 1)},
                     cv=5,
                     n_jobs=n_cpu-1,
                     verbose=3
                    )

# fcNN #
fcnn = fcNN(num_features=gene_expression_data.shape[1])
