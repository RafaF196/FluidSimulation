package collisions;

import org.lwjgl.util.vector.Vector3f;

public class TriangleCol {
	
	private Vector3f P1;
	private Vector3f P2;
	private Vector3f P3;
	private Vector3f normal;
	private Float d;
	
	public TriangleCol(Vector3f norm, Float d, Vector3f P1, Vector3f P2, Vector3f P3) {
		this.P1 = P1;
		this.P2 = P2;
		this.P3 = P3;
		this.normal = norm;
		this.d = d;
	}
	
	public Vector3f getP1() {
		return P1;
	}
	
	public Vector3f getP2() {
		return P2;
	}

	public Vector3f getP3() {
		return P3;
	}
	
	public Vector3f getNormal() {
		return normal;
	}
	
	public Float getD() {
		return d;
	}

}
