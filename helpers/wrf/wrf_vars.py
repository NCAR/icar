import operator as op

#  note, the elements of wrfvars must either be strings, operators, numbers, or lists.
#   lists will have their elements loaded one by one
#   strings will be treated as variable names to load from the file
#   operators will cause the next element to be loaded operated on with the last element
#   numbers will simply be used as is.
#
# Example : ["P",op.add,"PB"] means that the variable "P" will be loaded, then the variable "PB"
#           the result will be stored in variable "P" in the output file.


steps_per_day=24 # number of output steps per day in the WRF out file

var_conversions={"PH":"Z", "P":"pressure", "T":"potential_temperature"}

def rename_var(inputvar,dummy):
    """simple function to use as an operator to rename variables"""
    inputvar.name=var_conversions[inputvar.name]

wrfvars=["XLAT",   "XLONG",
         "XLAT_U", "XLONG_U",
         "XLAT_V", "XLONG_V",
         ["P",op.add,"PB", rename_var,0],
         [["PH",op.add,"PHB"], op.truediv, 9.8, rename_var, 0],
         "PSFC","HGT",
         "U","V",["T", op.add, 300, rename_var, 0],"QVAPOR",
         ["QCLOUD", op.add, "QRAIN"],
         ["QICE", op.add, "QSNOW", op.add, "QGRAUP"],
         "SWDOWN","GLW",
         "HFX","LH", "PBLH"
         ]

tsvar    = "TSK"
landmask = "XLAND"
