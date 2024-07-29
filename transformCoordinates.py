"""
Transform coordinate systems according to ICES SONAR-netcdf4 section 5.
To go from BCS to SCS, apply functions sequentially.
"""

import numpy as np

def transformBCStoTCS(posXYZ_BCS,anglePhiThetaPsi):
    # section 5.3.2, i.e. no stabilization is assumed.
    # No offsets between BCS and TCS
    # Typically the angles between BCS and TCS are all zero
    Phi = anglePhiThetaPsi[0]
    Theta = anglePhiThetaPsi[1]
    Psi = anglePhiThetaPsi[2]
    R_z_Psi = np.array([[np.cos(Psi), -np.sin(Psi), 0],
            [np.sin(Psi), np.cos(Psi),  0], 
            [0, 0, 1]
            ])
    R_y_Theta =np.array([[np.cos(Theta), 0, np.sin(Theta)],
            [0, 1, 0], 
            [-np.sin(Theta), 0, np.cos(Theta)]
            ])
    R_x_Phi =np.array([[1, 0, 0],
            [0, np.cos(Phi),  -np.sin(Phi)],
            [0, np.sin(Phi),  np.cos(Phi)]
            ])
    R_BCStoTCS = R_z_Psi.dot(R_y_Theta.dot(R_x_Phi))

    posXYZ_TCS = R_BCStoTCS.dot(posXYZ_BCS)
    return posXYZ_TCS


def transformTCStoPCS(posXYZ_TCS,angleXYZ,offsetXYZ,beamStabilized=False):
    # Compensate for angles and offsets between transducer face and platform coordinate system.
    # Current version assumes no stabilization of transducer beam.
    # posXYZ_TCS is the XYZ position in TCS coordinates
    # angleXYZ is the transducer tilt angles about the XYZ axis in PCS coordinates
    # offsetXYZ is the transducer offsets from the platform coordinates
    if beamStabilized:
        # implement stabilized algorithm
        print("Stabilized version pending implementation")
    elif not(beamStabilized):
        angleX = angleXYZ[0]
        angleY = angleXYZ[1]
        angleZ = angleXYZ[2]
        R_z_z = np.array([[np.cos(angleZ), -np.sin(angleZ), 0],
                [np.sin(angleZ), np.cos(angleZ),  0], 
                [0, 0, 1]
                ])
        R_y_y =np.array([[np.cos(angleY), 0, np.sin(angleY)],
                [0, 1, 0], 
                [-np.sin(angleY), 0, np.cos(angleY)]
                ])
        R_x_x =np.array([[1, 0, 0],
                [0, np.cos(angleX),  -np.sin(angleX)],
                [0, np.sin(angleX),  np.cos(angleX)]
                ])

        R_TCStoPCS = R_z_z.dot(R_y_y.dot(R_x_x))

        posXYZ_PCS = R_TCStoPCS.dot(posXYZ_TCS) + offsetXYZ
        return posXYZ_PCS

def transformPCStoSCS(posXYZ_PCS,angleHPR,offsetXYZ):
    # Compensates for platform heave, pitch and roll
    h = angleHPR[0]
    p = angleHPR[1]
    r = angleHPR[2]
    R_z_h = np.array([[np.cos(h), -np.sin(h), 0],
            [np.sin(h), np.cos(h),  0], 
            [0, 0, 1]
            ])
    R_y_p =np.array([[np.cos(p), 0, np.sin(p)],
            [0, 1, 0], 
            [-np.sin(p), 0, np.cos(p)]
            ])
    R_x_r =np.array([[1, 0, 0],
            [0, np.cos(r),  -np.sin(r)],
            [0, np.sin(r),  np.cos(r)]
            ])
    
    R_PCStoSCS = R_z_h.dot(R_y_p.dot(R_x_r))

    posXYZ_SCS = R_PCStoSCS.dot(posXYZ_PCS) + offsetXYZ
    return posXYZ_SCS

def transformSCStoTGS():
    # Transform from SCS coordinates  Terrestrial Geographical System
    print("Not yet implemented")

# Test:
angleBCS = [0,0,np.pi/2]
XYZ_BCS = [3,4,20]
XYZ_TCS = transformBCStoTCS(XYZ_BCS,angleBCS)
print("BCS coordinates: ",XYZ_BCS)
print("TCS coordinates: ",XYZ_TCS)

angleTCS = [-5*np.pi/180,0*5*np.pi/180,0]
offsetTCS = [1,-1,1]
XYZ_PCS = transformTCStoPCS(XYZ_TCS,angleTCS,offsetTCS)

print("PCS coordinates: ",XYZ_PCS)
