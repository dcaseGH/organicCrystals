import numpy as np

def angleBetweenVectors(v1, v2):
    return np.arccos(np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2)))

def numpyArrayToMathematicaString(npArray):
    ''' Copy this and paste into mathematica
        Assume matrix '''
    return '{{' + '},\n{'.join([','.join([str(c) for c in row]) for row in npArray]) + '}}'
        
def rmsdNoRestriction(pts1, pts2):
    ''' Find nearest neighbour for all points '''
    return (sum([min([np.linalg.norm(pts1[i] - x)**2 for x in pts2]) for i in xrange(pts1.shape[0])]) / pts1.shape[0])**0.5

def rmsd(pts1, pts2):
    ''' assume np arrays with shape npts,ndim '''
    assert(pts1.shape == pts2.shape)
    return (sum([np.linalg.norm(pts1[i] - pts2[i])**2 for i in xrange(pts1.shape[0])]) / pts1.shape[0])**0.5

def quaternion_axis_angle(axis, angle):
    #added 170417
    ax, ay, az = np.array(axis)/np.linalg.norm(axis)
    return np.array([     np.cos(angle/2),
                     ax * np.sin(angle/2),
                     ay * np.sin(angle/2),
                     az * np.sin(angle/2)])


def quaternion_mult ( q1, q2):
#    q3=range(0,4)
    q3 = np.zeros(4)
    q3[0] = q1[0] * q2[0] - q1[1] * q2[1] - q1[2] * q2[2] - q1[3] * q2[3]
    q3[1] = q1[0] * q2[1] + q1[1] * q2[0] + q1[2] * q2[3] - q1[3] * q2[2]
    q3[2] = q1[0] * q2[2] + q1[2] * q2[0] + q1[3] * q2[1] - q1[1] * q2[3]
    q3[3] = q1[0] * q2[3] + q1[3] * q2[0] + q1[1] * q2[2] - q1[2] * q2[1]
    return q3

def quaternion_rotatn(vec,quat):
#    dhc- this has a silly name,, but the old routine with a similar name is still used by pete
#    quat_inv = [quat[0],-1.0*quat[1],-1.0*quat[2],-1.0*quat[3]]
    quat_inv = np.array([quat[0],-1.0*quat[1],-1.0*quat[2],-1.0*quat[3]])
    return quaternion_mult( quaternion_mult(quat, np.array([0., vec[0], vec[1], vec[2]])), quat_inv )[1:4]

def list3norm(a=[None,None,None],b=[None,None,None]):
    import math
    diff=list3diff(a,b)
    return math.sqrt(math.pow(diff[0],2)+math.pow(diff[1],2)+math.pow(diff[2],2))

def list3normalize(a=[None,None,None]):
    n=list3norm(a,[0.0,0.0,0.0])
    return [a[0]/n,a[1]/n,a[2]/n]

def list3diff(a=[None,None,None],b=[None,None,None]):
    return [a[0]-b[0],a[1]-b[1],a[2]-b[2]]

def list3add(a=[None,None,None],b=[None,None,None]):
    return [a[0]+b[0],a[1]+b[1],a[2]+b[2]]

def quat_mult3 (q1, q2, q3):
    res=[0,0,0,0]
    t4=q1[0]
    t2=q2[0]
    t1=q1[1]
    t5=q2[1]
    t10=q1[2]
    t8=q2[2]
    t7=q1[3]
    t11=q2[3]
    t30=t10*t2
    t31=t5*t7
    t32=t4*t8
    t33=-t1*t11
    t34=t30+t31+t32+t33
    t14=q3[0]
    t23=t2*t7
    t24=-t10*t5
    t25=t1*t8
    t26=t11*t4
    t27=t23+t24+t25+t26
    t21=q3[1]
    t16=t2*t4
    t17=-t1*t5
    t18=-t10*t8
    t19=-t11*t7
    t20=t16+t17+t18+t19
    t28=q3[2]
    t3=t1*t2
    t6=t4*t5
    t9=-t7*t8
    t12=t10*t11
    t13=t12+t3+t6+t9
    t35=q3[3]
    res[1]=t13*t14+t20*t21-t27*t28+t34*t35
    res[2]=t21*t27+t20*t28+t14*t34-t13*t35
    res[3]=t14*t27+t13*t28-t21*t34+t20*t35
    return res

def rotation_by_q(vecin,theta,omega):
    from math import sin
    res = [0,0,0]
    try:
        res=[0,0,0]
        t1=omega[0]
        t2=t1*t1
        t3=omega[1]
        t4=t3*t3
        t5=omega[2]
        t6=t5*t5
        t7=t2+t4+t6
        t12=5.e-1*theta
        t13=sin(t12)
        t14=t13*t13
        t15=1.0/t14
        t9=4.*t14
        t10=sin(theta)
        t16=vecin[0]
        t19=vecin[1]
        t21=vecin[2]
        t8=2.5e-1/t7
        t18=t7**5.e-1
        t30=t19*t3
        res[0]=(t10*(t10*t7*t15*t16-4.*t18*(-t21*t3+t19*t5))+t9*(t16*t2+2.e0*t1*(t30+t21*t5)-t16*(t4+t6)))*t8
        res[1]=(t9*(-t19*t2+2.e0*t1*t16*t3+t19*(t4-t6)+2.*t21*t3*t5)+t10*(t10*t7*t15*t19+4.e0*t18*(-t1*t21+t16*t5)))*t8
        res[2]=(t10*(t10*t7*t15*t21-4.*t18*(-t1*t19+t16*t3))+t9*(-t21*(t2+t4-t6)+2.e0*(t1*t16+t30)*t5))*t8
    except ZeroDivisionError:
        return vecin
    return res

def rotation_by_q_old(vec_in,theta,omega):
    import math

    rvnorm = (omega[0]**2 + omega[1]**2 + omega[2]**2)**0.5

    st2 = math.sin(theta/2.0) / rvnorm
    quat_rot = [math.cos( theta/2.0 ) , omega[0] * st2, omega[1] * st2, omega[2] * st2]
    quat_rot_inv = [quat_rot[0],-quat_rot[1],-quat_rot[2],-quat_rot[3]]
    temp_quat2 = [0.0 , vec_in[0] , vec_in[1] , vec_in[2]]
    vec_out = quaternion_mult( quaternion_mult(quat_rot, temp_quat2), quat_rot_inv )

    return vec_out[1:4]

def list3cross(a=[None,None,None],b=[None,None,None]):
    c = [a[1]*b[2] - a[2]*b[1],
         a[2]*b[0] - a[0]*b[2],
         a[0]*b[1] - a[1]*b[0]
         ]
    return c

def list3dot(a=[None,None,None],b=[None,None,None]):
    return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]

def list3mag(a=[None,None,None]):
    import math
    return math.sqrt(math.pow(a[0],2)+math.pow(a[1],2)+math.pow(a[2],2))


def list3subtract(a=[None,None,None],b=[None,None,None]):
    c = (a[0]-b[0],a[1]-b[1],a[2]-b[2])
    return c

def list3multiply(a=[None,None,None],fac=None):
    return (fac*a[0],fac*a[1],fac*a[2])

#    def list3vecmatrix(v=[None,None,None],m=[[None,None,None],[None,None,None],[None,None,None]]):
#        return [v[0]*m[0][0]+v[0]*m[0][0],0,0]

def list3vecmatrix(a=[None,None,None],m=[[None,None,None],[None,None,None],[None,None,None]]):
    v1=list3multiply(m[0],a[0])
    v2=list3multiply(m[1],a[1])
    v3=list3multiply(m[2],a[2])
    return list3add(list3add(v1,v2),v3)

def list3_proj_pts_axis ( axis , set_vertices , set_radii):
    proj = list3dot(axis , set_vertices[0]  )
    my_min = proj - set_radii[0]
    my_max = proj + set_radii[0]
    for i in xrange ( 1, len(set_vertices) ) :
        proj = list3dot(axis , set_vertices[i] )
        if proj - set_radii [ i ] < my_min :
            my_min = proj - set_radii [ i ]
        if proj + set_radii [ i ] > my_max :
            my_max = proj + set_radii [ i ]
    return my_max , my_min

def list3delparallelvecs(set_vecs , param):
    len_vecs = len(set_vecs)
    i,j=0,0
    while i < len_vecs :
        j=i+1
        while j < len_vecs :
#                if abs(list3dot (set_vecs[i],set_vecs[j])/(list3mag(set_vecs[i])*list3mag(set_vecs[j]))) > param:
            if abs(list3dot (set_vecs[i],set_vecs[j])) > param:
                del set_vecs[j]
                len_vecs-=1
            else:
                j+=1
        i+=1
    return set_vecs

def list3matrixmatrix(m=[[None,None,None],[None,None,None],[None,None,None]],n=[[None,None,None],[None,None,None],[None,None,None]]):
    return (
            (m[0][0]*n[0][0]+m[0][1]*n[1][0]+m[0][2]*n[2][0],
             m[0][0]*n[0][1]+m[0][1]*n[1][1]+m[0][2]*n[2][1],
             m[0][0]*n[0][2]+m[0][1]*n[1][2]+m[0][2]*n[2][2]),
            (m[1][0]*n[0][0]+m[1][1]*n[1][0]+m[1][2]*n[2][0],
             m[1][0]*n[0][1]+m[1][1]*n[1][1]+m[1][2]*n[2][1],
             m[1][0]*n[0][2]+m[1][1]*n[1][2]+m[1][2]*n[2][2]),
            (m[2][0]*n[0][0]+m[2][1]*n[1][0]+m[2][2]*n[2][0],
             m[2][0]*n[0][1]+m[2][1]*n[1][1]+m[2][2]*n[2][1],
             m[2][0]*n[0][2]+m[2][1]*n[1][2]+m[2][2]*n[2][2]))

def list3make_all_translations ( lat_vecs , depth ):
    vector_list = []
    for ia in range (-1 * depth , depth + 1 ):
        for ib in range (-1 * depth , depth + 1 ):
            for ic in range (-1 * depth , depth + 1 ):
                 vector_list.append(list3vecmatrix([ia,ib,ic] , lat_vecs))
    return sorted(vector_list, key = lambda x: x[0]**2 + x[1]**2 + x[2]**2 )

def list3cross_combinations(set_vecs1 , set_vecs2):
    set_out = []
    for i in range ( 0 , len(set_vecs1) ):
        for j in range ( 0 , len(set_vecs2) ):
            try:
                if abs(list3dot(list3normalize(set_vecs1[i]),list3normalize(set_vecs2[j]))) < 1.0:
                    set_out.append(list3normalize(list3cross( set_vecs1[i] , set_vecs2[j] )))
            except:
                pass
    return set_out

def list3cross_comb_tri(set_vecs1 , set_vecs2):
    set_out = []
    for i in range ( 0 , len(set_vecs1) ):
        for j in range ( 0 , i ):
            try:
                if abs(list3dot(list3normalize(set_vecs1[i]),list3normalize(set_vecs2[j]))) < 1.0:
                    set_out.append(list3normalize(list3cross( set_vecs1[i] , set_vecs2[j] )))
            except:
                pass
    return set_out

def list3edge_vecs(edges,op_mat,op_vec):
    vecs_edges = []
    for i in range ( 0 , len(edges) ) :
        for j in range ( 0 , i ) : #NOTE that for op= x,y,z we may have trouble when going back to above line as edges overlap
            vecs_edges.append(list3normalize( list3cross( list3subtract( edges[i][1],edges[i][0] ),list3subtract(list3add(list3vecmatrix(edges[j][1],op_mat),op_vec),list3add(list3vecmatrix(edges[j][0],op_mat),op_vec)))))
    return vecs_edges

def list3perp_face_vecs(faces):
    perp_vectors = []
    for face_pt in faces:
        perp_vectors.append(list3normalize ( list3cross(list3subtract(face_pt[2],face_pt[1]),list3subtract(face_pt[1],face_pt[0]))))
    return perp_vectors

def list3list_add(list_vecs, add_vec):
    return [list3add(vec,add_vec)  for vec in list_vecs]

def list3apply_sym(vec , op_mat , op_vec):
    return list3add(list3vecmatrix(vec,op_mat),op_vec)

def list3apply_sym_list_norm(vec_list , op_mat , op_vec):
    new_vec_list = []
    for vec in vec_list:
        new_vec_list.append(list3normalize(list3add(list3vecmatrix(vec,op_mat),op_vec)))
    return new_vec_list

def list3apply_sym_list(vec_list , op_mat , op_vec):
    new_vec_list = []
    for vec in vec_list:
        new_vec_list.append(list3add(list3vecmatrix(vec,op_mat),op_vec))
    return new_vec_list

def list3required_dvec( overlap , vector , axis ):
    cos_rdv_theta =  list3dot(vector , axis) / ( list3mag(vector) * list3mag (axis) )
    disp_vec = list3multiply ( list3normalize(axis) , (overlap / (cos_rdv_theta + 0.0000001)) )
    return disp_vec

def list3find_overlap( EP1, EP2):
    EP10, EP11 = EP1
    EP20, EP21 = EP2
    if( EP10 > EP20 ):
        if( EP11 < EP20):
            if( EP11 > EP21):
                return abs(EP20 - EP11)
            else:
                if abs(EP10 - EP21) <= abs(EP11 -EP20):
                    return abs(EP10 - EP21)
                else:
                    return abs(EP11 - EP20)
        else:
            return 0.0
    else:
        if( EP10 > EP21 ):
            if( EP11 < EP21):
                return abs(EP10 - EP21)
            else:
                if abs(EP10 - EP21) <= abs(EP11 -EP20):
                    return abs(EP10 - EP21)
                else:
                    return abs(EP11 - EP20)
        else:
             return 0.0

def list3list_dot( l1, l2):
    import numpy as np
    return list3list_dot(np.array(l1), np.array(l2))

def list3square_trace(l1):
    # Tr (l1 ^T . l1)
    return sum([list3dot(x, x) for x in l1])

def optimalZRotation(pts1, pts2, allowRMSDNoRestriction = False):
    ''' Find best overlay (1 onto 2 -and give quaternion for this rotation) for angle around z '''
    from scipy.optimize import fmin #clearly analytic result must be possible #powells method may be better?                                                                         
    from quaternions import quaternion_rotatn, make_quat_rot2
#    from listMathsSlow import rmsd

    def functionOfOverlay(angle, pts1, pts2, allowRMSDNoRestriction):
        _quat = make_quat_rot2(angle, np.array([0., 0., 1.]))
        rotatedPts1 = np.array([quaternion_rotatn(x, _quat) for x in pts1])
        if allowRMSDNoRestriction:
            return rmsdNoRestriction(rotatedPts1, pts2)
        else:
            return rmsd(rotatedPts1, pts2)

    optimisationOutputs = fmin(functionOfOverlay, 0., (pts1, pts2, allowRMSDNoRestriction), full_output=True, disp=False)

    return (optimisationOutputs[1],
            make_quat_rot2(optimisationOutputs[0][0], np.array([0., 0., 1.])))
            
#    return (make_quat_rot2(optimisationOutputs[0][0], np.array([0., 0., 1.])),
#            optimisationOutputs[1])


def overlay_points_RMSD(set_pts1, set_pts2):
    ''' This returns the optimal RMSD, and the quaternion to bring about rotation of set 2 onto set 1
    Assumes that both sets have mean point at origin
    See DOI 10.1002/jcc.20110 for theory and naming conventions
     '''

    #N.B. ORIGINALLY BELIEVED TO ROTATE 2 ONTO 1 BUT THIS DOES NOT SEE TO BE THE CASE! OTHER WAY AROUND (SEE TESTS)
    import numpy as np
    set_points1 = np.array(set_pts1)
    set_points2 = np.array(set_pts2)
    if False:
        matrix_R = np.array(list3list_dot(set_points1, set_points2))
    else:
        matrix_R = np.dot(np.transpose(set_points1), set_points2)
#    print matrix_R, np.dot(np.transpose(set_points1), set_points2)
    matrix_F = np.array([[matrix_R[0,0] + matrix_R[1,1] + matrix_R[2,2],\
                             matrix_R[1,2] - matrix_R[2,1],\
                             matrix_R[2,0] - matrix_R[0,2],\
                             matrix_R[0,1] - matrix_R[1,0]],\
                         [matrix_R[1,2] -matrix_R[2,1],\
                             matrix_R[0,0] - matrix_R[1,1] - matrix_R[2,2],\
                             matrix_R[0,1] + matrix_R[1,0],\
                             matrix_R[0,2] + matrix_R[2,0]],\
                         [matrix_R[2,0] - matrix_R[0,2],\
                             matrix_R[0,1] + matrix_R[1,0],\
                             -matrix_R[0,0] + matrix_R[1,1] - matrix_R[2,2],\
                             matrix_R[1,2] + matrix_R[2,1]],\
                         [matrix_R[0,1] - matrix_R[1,0],\
                             matrix_R[0,2] + matrix_R[2,0],\
                             matrix_R[1,2] + matrix_R[2,1],\
                             -matrix_R[0,0] - matrix_R[1,1] + matrix_R[2,2]]])
    eigenValues,eigenVectors = np.linalg.eig(matrix_F)
    eigenVectors = [x for (y,x) in sorted(zip(eigenValues , eigenVectors.transpose()))]
    eigenValues  = sorted(eigenValues)
    def best_fit_MSD(l1, l2, max_lambda):
#        return (list3square_trace(l1) + list3square_trace(l2) - 2.0 * max_lambda) / l1.shape[0]
        return ((sum([list3dot(x, x) for x in l1]) + sum([list3dot(x, x) for x in l2]) - 2.0 * max_lambda) / l1.shape[0])
    MSD = best_fit_MSD(set_points1, set_points2, eigenValues[-1])
    if MSD > 0.0:
        return (MSD ** 0.5, tuple(eigenVectors[-1]))
    elif MSD < -0.1:
        print "MSD is negative, %s, in overlay of points" %(MSD)
        return (MSD, tuple(eigenVectors[-1]))
    else:
        return (0.0, tuple(eigenVectors[-1]))

def map_into_hypersphere(v1):#r=1
    ''' maps n dimensional v1 element of [0,1)^n into Real^n, with Norm < r=1
        this method should be attributed to Roger Stafford, but cannot find a paper,
        only seen on Matlab forums. Naming conventions from here -dhc 190115 '''
    from scipy.stats import norm
    from scipy.special import gamma, gammainc
    n = float(len(v1))
    #norm.ppf is the percentile point function, i.e. generates normal distribution with mean = 0, sigma = 1
    v1 = map(norm.ppf, v1)
    s2 = (sum([x * x for x in v1]))
    inv_sqrt_s2 = s2 ** -0.5
    #gamma_inc is the (regularized??) incomplete lower Gamma function- don't divide by gamma (online scipy information is ambiguous, but this way works)
    scale = inv_sqrt_s2 * (gammainc(n / 2., s2 / 2.)) ** ( 1. / n )
    return [x * scale for x in v1]
