''' quaternions have order [scalar, vector[0], v[1], v[2]] '''

import numpy as np

def make_quat_rot1( x3):
    import math
    pi = math.acos(-1.0)
    x0 = x3[0];x1 = x3[1] ; x2 = x3[2]
    quat_rot = [(1.0 - x0)**0.5 * math.sin(2.0*pi*x1) ,(1.0 - x0)**0.5 * math.cos(2.0*pi*x1) ,x0**0.5 * math.sin(2.0*pi*x2) , x0**0.5 * math.cos(2.0*pi*x2)]
    return np.array(quat_rot)

def make_quat_rot2(theta , omega):
    import math
    quat_rot = [math.cos( theta/2.0 ) , omega[0] * math.sin(theta/2.0), omega[1] * math.sin(theta/2.0), omega[2] * math.sin(theta/2.0)]
    return np.array(quat_rot)

def inverse_quat_rot2(q):
    ''' quaternion -> angle, (normalized) axis '''
    _theta = np.arccos(q[0]) * 2.
    return _theta, q[1:]/ np.sin(_theta/2.)

def inverse_quat_rot2Safe(q):
    ''' quaternion -> angle, (normalized) axis '''
    _theta = np.arccos(q[0]) * 2.
    #set convention that 0. radian rotation is around z
    if _theta < 1.e-9:
        return 0., np.array([0., 0., 1.])
    return _theta, q[1:]/ np.sin(_theta/2.)

def zxz_quaternion(p, t, P):
    ''' (phi, theta, Phi) around z x z '''
    return np.array([np.cos((p+P)/2.) * np.cos(t/2.),
                     np.cos((p-P)/2.) * np.sin(t/2.),
                     np.sin((p-P)/2.) * np.sin(t/2.),
                     np.sin((p+P)/2.) * np.cos(t/2.)])

def zyz_quaternion(p, t, P):
    ''' (phi, theta, Phi) around z x z '''
    return np.array([np.cos((p+P)/2.) * np.cos(t/2.),
                    -np.sin((p-P)/2.) * np.sin(t/2.),
                     np.cos((p-P)/2.) * np.sin(t/2.),
                     np.sin((p+P)/2.) * np.cos(t/2.)])

def zxz_angles(q):
    ''' extrinsic 3-1-3 Euler angles 
    https://en.wikipedia.org/wiki/Rotation_formalisms_in_three_dimensions#Euler_angles_.28z-x-z_extrinsic.29_.E2.86.92_Quaternion 
    quaternion -> (p, t, P) '''

    return np.array([np.arctan2(q[1]*q[3] + q[2]*q[0], -1.*(q[2]*q[3] - q[1]*q[0])),
                     np.arccos(-q[1]**2 - q[2]**2 + q[3]**2 + q[0]**2),
                     np.arctan2(q[1]*q[3] - q[2]*q[0], q[2]*q[3] + q[1]*q[0])])

def zyz_angles(q):
    ''' found this by hacking '''
    return np.array([np.arctan2(q[2]*q[3] - q[1]*q[0], q[1]*q[3] + q[2]*q[0]),
                     np.arccos(-q[1]**2 - q[2]**2 + q[3]**2 + q[0]**2),
                     np.arctan2(q[2]*q[3] + q[1]*q[0], -1.*(q[1]*q[3] - q[2]*q[0]))])


def quatToMatrix(quat):
    ''' Unit quaternion to (left handed i.e. post multiplied) matrix
        i.e. matrix.vector = rotatedVector
    https://en.wikipedia.org/wiki/Rotation_matrix#General_rotations '''
    w, x, y, z = quat[0], quat[1], quat[2], quat[3]
    return np.array([[1. - 2.*y**2 - 2.*z**2, 2.*x*y - 2.*z*w,         2.*x*z + 2.*y*w      ],
                     [2.*x*y + 2.*z*w,        1. - 2.*x**2 - 2.*z**2,  2.*y*z - 2.*x*w      ],
                     [2.*x*z - 2.*y*w,        2.*y*z + 2.*x*w,         1. -2.*x**2 - 2.*y**2]])
                     
def quaternion_mult ( q1, q2):
    q3 = np.zeros(4)#range(0,4)
    q3[0] = q1[0] * q2[0] - q1[1] * q2[1] - q1[2] * q2[2] - q1[3] * q2[3]
    q3[1] = q1[0] * q2[1] + q1[1] * q2[0] + q1[2] * q2[3] - q1[3] * q2[2]
    q3[2] = q1[0] * q2[2] + q1[2] * q2[0] + q1[3] * q2[1] - q1[1] * q2[3]
    q3[3] = q1[0] * q2[3] + q1[3] * q2[0] + q1[1] * q2[2] - q1[2] * q2[1]
    return np.array(q3)

def quaternion_rotatn(vec,quat):
#    dhc- this has a silly name,, but the old routine with a similar name is still used by pete
    quat_inv = np.array([quat[0],-1.0*quat[1],-1.0*quat[2],-1.0*quat[3]])
    return quaternion_mult( quaternion_mult(quat, np.array([0.0, vec[0], vec[1], vec[2]])), quat_inv )[1:4]
