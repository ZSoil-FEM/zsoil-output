#binary codes for plastic behavior type
#combined plastic states are possible

class PlasticCodes ():

    PL_CODE_VIRG_EL         = 0
    PL_CODE_UNL_FROM_PL     = 1
    PL_CODE_SHEAR_CURR      = 2
    PL_CODE_SIG12           = 4
    PL_CODE_SIG23           = 8
    PL_CODE_SHEAR_LIM       = 16
    PL_CODE_CAP             = 32
    PL_CODE_TENS            = 256
    PL_CODE_TENS_CURR       = 512
    PL_CODE_SHEAR_PHYS_1    = 1024
    PL_CODE_SHEAR_PHYS_2    = 2048
    PL_CODE_SHEAR_PHYS_3    = 4096
    PL_CODE_TENS_PHYS_1     = 8192
    PL_CODE_TENS_PHYS_2     = 16384
    PL_CODE_TENS_PHYS_3     = 32768
    PL_CODE_TENS_POST_PEAK  = 65536
    PL_CODE_COMPR_POST_PEAK = 131072
    PL_CODE_SHEAR_POST_PEAK = 262144
    PL_CODE_COMPR           = 524288
    PL_CODE_COMPR_CURR      = 1048576

