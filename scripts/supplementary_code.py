import numpy as np
from sklearn.linear_model import LinearRegression

# logit regression definition
class LogitRegression(LinearRegression):

    def fit(self, x, p):
        p = np.asarray(p)
        y = np.log(p / (1 - p))
        return super().fit(x, y)

    def predict(self, x):
        y = super().predict(x)
        return 1 / (np.exp(-y) + 1)

# ctype cv generator
def ctype_cv_generator(gene_expr_train, types_data):
    
    gene_expr_train_t = gene_expr_train.join(types_data)
    types = set(gene_expr_train_t.type) # only go through 20 types that are in the training set
    
    for ctype in types:
        
        test_idx_names = gene_expr_train_t.loc[gene_expr_train_t.type == ctype, :].index
        test_idx_ints = gene_expr_train_t.index.get_indexer(test_idx_names)

        train_idx_names = gene_expr_train_t.index.drop(test_idx_names)
        train_idx_ints = gene_expr_train_t.index.get_indexer(train_idx_names)
        
        yield train_idx_ints, test_idx_ints

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
#                     kernel_regularizer=regularizers.l2(0.001)
                )
            )

#         # hidden layers
#         for i in range(1):
#             model.add(
#                 Dense(
#                     units=(units_init // 2),
#                     activation='relu',
#                     kernel_regularizer=regularizers.l2(0.01)
#                 ))

# #             model.add(Dropout(rate=0.1))

#         model.add(
#             Dense(
#                 units=(units_init // 4),
#                 activation='relu',
#                 kernel_regularizer=regularizers.l2(0.01)
#             ))

#         model.add(
#             Dense(
#                 units=(units_init // 8),
#                 activation='relu',
#                 kernel_regularizer=regularizers.l2(0.01)
#             ))

#         model.add(
#             Dense(
#                 units=(units_init // 16),
#                 activation='relu'))

        # output layer
        model.add(Dense(1, activation='sigmoid'))

        model.compile(
            optimizer=Adam(learning_rate/1000),
            loss='mean_squared_error')
        
        return model
    
    
### models definitions ###

# Linear Regression #
linreg = LinearRegression()

# Logit Reg #
model = LogitRegression()

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
                            'C':[0.001, 0.01, 0.1, 1, 10]}, # adjusted
                 cv=5,
                 factor=3,
                 scoring='neg_mean_squared_error', # changed
                 n_jobs=n_cpu-1,
                 verbose=3,
                 error_score='raise'
                )

# Gradient Boosting Regressor #
gradboost = HalvingGridSearchCV(estimator=GradientBoostingRegressor(), 
                    param_grid={'max_depth':np.arange(1, 101, 1)},
                     cv=5,
                     n_jobs=n_cpu-1,
                     verbose=3
                    )

# fcNN #
fcnn = fcNN(num_features=gene_expression_data.shape[1]) # epochs 6, batch size 1
