package toolbox;

import org.lwjgl.util.vector.Matrix4f;
import org.lwjgl.util.vector.Vector3f;
import org.lwjgl.util.vector.Vector4f;

import entities.Camera;

public class MatrixOps {

	public static Matrix4f createTransformationMatrix(Vector3f translation, float rx, float ry,
			float rz, float scale) {
		Matrix4f matrix = new Matrix4f();
		matrix.setIdentity();
		Matrix4f.translate(translation, matrix, matrix);
		Matrix4f.rotate((float) Math.toRadians(rx), new Vector3f(1,0,0), matrix, matrix);
		Matrix4f.rotate((float) Math.toRadians(ry), new Vector3f(0,1,0), matrix, matrix);
		Matrix4f.rotate((float) Math.toRadians(rz), new Vector3f(0,0,1), matrix, matrix);
		Matrix4f.scale(new Vector3f(scale,scale,scale), matrix, matrix);
		return matrix;
	}
	
	public static Matrix4f createViewMatrix(Camera camera) {
		Matrix4f viewMatrix = new Matrix4f();
		viewMatrix.setIdentity();
		Matrix4f.rotate((float) Math.toRadians(camera.getPitch()), new Vector3f(1, 0, 0), viewMatrix,
				viewMatrix);
		Matrix4f.rotate((float) Math.toRadians(camera.getYaw()), new Vector3f(0, 1, 0), viewMatrix, viewMatrix);
		Vector3f cameraPos = camera.getPosition();
		Vector3f negativeCameraPos = new Vector3f(-cameraPos.x,-cameraPos.y,-cameraPos.z);
		Matrix4f.translate(negativeCameraPos, viewMatrix, viewMatrix);
		return viewMatrix;
	}
	
	// Compute Qv for a vertex (method used on the slides)
	public static Matrix4f createQuadricMatrix(Vector3f pos, Vector3f norm) {
		Matrix4f qvMat = new Matrix4f();
		Vector4f pv = new Vector4f(norm.x, norm.y, norm.z, -Vector3f.dot(norm, pos));
		qvMat.m00 = pv.x*pv.x; qvMat.m10 = pv.x*pv.y; qvMat.m20 = pv.x*pv.z; qvMat.m30 = pv.x*pv.w;
		qvMat.m01 = pv.y*pv.x; qvMat.m11 = pv.y*pv.y; qvMat.m21 = pv.y*pv.z; qvMat.m31 = pv.y*pv.w;
		qvMat.m02 = pv.z*pv.x; qvMat.m12 = pv.z*pv.y; qvMat.m22 = pv.z*pv.z; qvMat.m32 = pv.z*pv.w;
		qvMat.m03 = pv.w*pv.x; qvMat.m13 = pv.w*pv.y; qvMat.m23 = pv.w*pv.z; qvMat.m33 = pv.w*pv.w;
		return qvMat;
	}
	
	// Compute the optimal position of a vertex given its Q matrix
	public static Vector3f obtainOptimalContractionVertex(Matrix4f qMat) {
		Vector3f newVer = new Vector3f();
		
		qMat.m03 = 0; qMat.m13 = 0; qMat.m23 = 0; qMat.m33 = 1;
		// Only invert if the determinant is not 0 (non-invertible otherwise)
		if (qMat.determinant() != 0) {
			qMat.invert();
			// Results of multiplying with column vector [0,0,0,1]
			newVer.x = qMat.m30;
			newVer.y = qMat.m31;
			newVer.z = qMat.m32;
		} else { // If we can't invert the matrix, try another method
			return null;
		}
		
		return newVer;
	}
	
}
