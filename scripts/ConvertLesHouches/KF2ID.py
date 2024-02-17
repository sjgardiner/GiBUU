"""
KF to ID scheme

cf. ID_translation/kfbuuC:

"""


KF_CODES = {
    # non strange baryons:
    2112: [  1, 0], # p,n
    2212: [  1, 1],
    1114: [  2,-1], # Delta
    2114: [  2, 0],
    2214: [  2, 1],
    2224: [  2, 2],
    # strange baryons:
    3122: [ 32, 0], # Lambda
    3112: [ 33,-1], # Sigma
    3212: [ 33, 0],
    3222: [ 33, 1],
    3114: [ 34,-1], # Sigma*
    3214: [ 34, 0],
    3224: [ 34, 1],
    # doule strange baryons:
    3312: [ 53,-1], # Xi
    3322: [ 53, 0],
    3314: [ 54,-1], # Xi*
    3324: [ 54, 0],
    # triple strange baryons:
    3334: [ 55,-1], # Omega


    # non strange mesons:
    -211: [101, -1], # pi
     111: [101,  0],
     211: [101,  1],
     221: [102,  0], # eta
    -213: [103, -1], # rho
     113: [103,  0],
     213: [103,  1],
     223: [105,  0], # omega
     331: [106,  0], # eta prime
     333: [107,  0], # phi
    # strange mesons:
     311: [110,  1], # K+
     321: [110,  0], # K0
    -321: [111, -1], # K-
    -311: [111,  0], # K0bar

     313: [112,  1], # K+*
     323: [112,  0], # K0*
    -323: [113, -1], # K-*
    -313: [113,  0], # K0bar*

    }


# Complete the missing antiparticles
for i in list(KF_CODES.keys()):
    if i > 0 and -i not in KF_CODES:
        KF_CODES[-i] = [-KF_CODES[i][0],-KF_CODES[i][1]]

del i
