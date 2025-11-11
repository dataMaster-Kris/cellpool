#importing cellpose downloads more files from the web
#Compute clusters may not provide internet access on the compute nodes
#Hence, the first time import of cellpose should happen during container creation
from cellpose import core, utils, io, models, metrics
import numpy as np

mdl = models.CellposeModel(model_type='tissuenet')

a = np.ones((500, 500))
mdl.eval(a) 
