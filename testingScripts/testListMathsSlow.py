import numpy as np
import unittest
import os

class TestOverlays(unittest.TestCase):

    def setUp(self):
        pass

    def test_randomStuff(self):
        from listMathsSlow import angleBetweenVectors
        self.assertEqual(angleBetweenVectors(np.array([0., 0., 1.]),
                                             np.array([0., 0., -1000.])),
                         np.pi)
        self.assertEqual(angleBetweenVectors(np.array([0., 0., 1.]),
                                             np.array([0., -1000., 0.])),
                         np.pi/2.)

    def test_optimisingOverlays(self):

        pts1 = np.array([[-0.541357990208,0.840223963793,0.967951435],
			 [-0.507301547591,-0.542664198999,1.171816351],
			 [0.0255987444643,-1.38443421464,0.204678748],
			 [0.556845352422,-0.849058545788,-0.964797836],
			 [0.509212767371,0.526541086859,-1.159201955],
			 [-0.0429973264567,1.40939190877,-0.220446743]])

	pts2 = np.array([[-0.452870161944,1.0471304112,0.822046103134],
		         [0.838143143533,0.525685681579,0.978663279603],
			 [1.28692322978,-0.496230348673,0.148563133317],
			 [0.457098160295,-1.03724996618,-0.829129779641],
			 [-0.837847702903,-0.538192702582,-0.964350371592],
	        	 [-1.29144666876,0.49885692466,-0.155792364821]])

        np.testing.assert_array_almost_equal(np.mean(pts1, axis=0),
                                             np.mean(pts2, axis=0))

	from listMathsSlow import rmsd
	print rmsd(pts1, pts2)

	from listMathsSlow import overlay_points_RMSD
        originalRMSD, qRot = overlay_points_RMSD(pts1, pts2)
#	print originalRMSD
        from quaternions import quaternion_rotatn
        self.assertAlmostEqual(rmsd(pts2,
                                    np.array([quaternion_rotatn(x, qRot) for x in pts1])),
                               originalRMSD)

        from listMathsSlow import optimalZRotation
        bestRMSDZ, bestQRotZ = optimalZRotation(pts1,
                                                pts2)

        from quaternions import make_quat_rot2
        from listMathsSlow import rmsdNoRestriction

        #check with a grid that optimal Z really is about the best
        _nPts = 20
        gridOptimisedPoints = np.array([rmsd(pts2,
                                             np.array([quaternion_rotatn(x, 
                                                                         make_quat_rot2(float(i)/float(_nPts) * 2. * np.pi,
                                                                                        np.array([0., 0., 1.]))) 
                                                       for x in pts1])
                                             )
                                        for i in xrange(_nPts)])
        self.assertTrue(all([x > bestRMSDZ for x in gridOptimisedPoints]))
        self.assertTrue(bestRMSDZ + 1.e-3 > gridOptimisedPoints[3])



if __name__ == '__main__':
    unittest.main()
