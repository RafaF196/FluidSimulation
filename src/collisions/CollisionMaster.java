package collisions;

import java.util.ArrayList;

import org.lwjgl.util.vector.Vector3f;

import models.RawModel;
import objFileLoader.ObjFileLoader;
import render.Loader;

public class CollisionMaster {
	
	private ArrayList<PlaneCol> planes;
	private ArrayList<TriangleCol> triangles;
	private ArrayList<SphereCol> spheres;
	
	public CollisionMaster () {
		this.planes = new ArrayList<PlaneCol>();
		this.triangles = new ArrayList<TriangleCol>();
		this.spheres = new ArrayList<SphereCol>();
	}
	
	// Add colision info and generate the model
	public RawModel addPlane(Loader loader, Vector3f ori, Float val) {
		float[] vArray = new float[12];
		float[] tArray = new float[8];
		float[] nArray = new float[12];
		int[] iArray = new int[6];
		vArray[0] = -40; vArray[1] = -val; vArray[2] = -40;
		vArray[3] = 40; vArray[4] = -val; vArray[5] = -40;
		vArray[6] = -40; vArray[7] = -val; vArray[8] = 40;
		vArray[9] = 40; vArray[10] = -val; vArray[11] = 40;
		tArray[0] = 0; tArray[1] = 0; tArray[2] = 0; tArray[3] = 1;
		tArray[4] = 1; tArray[5] = 0; tArray[6] = 1; tArray[7] = 1;
		nArray[0] = 0; nArray[1] = 1; nArray[2] = 0;
		nArray[3] = 0; nArray[4] = 1; nArray[5] = 0;
		nArray[6] = 0; nArray[7] = 1; nArray[8] = 0;
		nArray[9] = 0; nArray[10] = 1; nArray[11] = 0;
		iArray[0] = 1; iArray[1] = 2; iArray[2] = 3;
		iArray[3] = 2; iArray[4] = 1; iArray[5] = 0;
		planes.add(new PlaneCol(ori, val));
		return loader.loadToVAO("plane", vArray, tArray, nArray, iArray);
	}
	
	public RawModel addTriangle(Loader loader, Vector3f P1, Vector3f P2, Vector3f P3) {
		Vector3f V = new Vector3f();
		Vector3f U = new Vector3f();
		Vector3f ori = new Vector3f();
		Vector3f.sub(P2, P1, U);
		Vector3f.sub(P3, P1, V);
		Vector3f.cross(U, V, ori);
		ori.normalise();
		float[] vArray = new float[9];
		float[] tArray = new float[6];
		float[] nArray = new float[9];
		int[] iArray = new int[6];
		vArray[0] = P1.x; vArray[1] = P1.y; vArray[2] = P1.z;
		vArray[3] = P2.x; vArray[4] = P2.y; vArray[5] = P2.z;
		vArray[6] = P3.x; vArray[7] = P3.y; vArray[8] = P3.z;
		tArray[0] = 0; tArray[1] = 0; tArray[2] = 0; tArray[3] = 1;
		tArray[4] = 1; tArray[5] = 0;
		nArray[0] = ori.x; nArray[1] = ori.y; nArray[2] = ori.z;
		nArray[3] = ori.x; nArray[4] = ori.y; nArray[5] = ori.z;
		nArray[6] = ori.x; nArray[7] = ori.y; nArray[8] = ori.z;
		iArray[0] = 0; iArray[1] = 1; iArray[2] = 2;
		iArray[3] = 1; iArray[4] = 0; iArray[5] = 2;
		triangles.add(new TriangleCol(ori, -Vector3f.dot(ori, P1), P1, P2, P3));
		return loader.loadToVAO("triangle", vArray, tArray, nArray, iArray);
	}
	
	public RawModel addSphere(Loader loader, Vector3f center, Float radius) {
		spheres.add(new SphereCol(center, radius));
		return ObjFileLoader.loadOBJ("uvsphere", loader);
	}
	
	// Used for models
	public void addTriangleNoLoad(Vector3f P1, Vector3f P2, Vector3f P3) {
		Vector3f V = new Vector3f();
		Vector3f U = new Vector3f();
		Vector3f ori = new Vector3f();
		Vector3f.sub(P2, P1, U);
		Vector3f.sub(P3, P1, V);
		Vector3f.cross(U, V, ori);
		if (ori.length() == 0) return;
		ori.normalise();
		float[] vArray = new float[9];
		float[] tArray = new float[6];
		float[] nArray = new float[9];
		int[] iArray = new int[6];
		vArray[0] = P1.x; vArray[1] = P1.y; vArray[2] = P1.z;
		vArray[3] = P2.x; vArray[4] = P2.y; vArray[5] = P2.z;
		vArray[6] = P3.x; vArray[7] = P3.y; vArray[8] = P3.z;
		tArray[0] = 0; tArray[1] = 0; tArray[2] = 0; tArray[3] = 1;
		tArray[4] = 1; tArray[5] = 0;
		nArray[0] = ori.x; nArray[1] = ori.y; nArray[2] = ori.z;
		nArray[3] = ori.x; nArray[4] = ori.y; nArray[5] = ori.z;
		nArray[6] = ori.x; nArray[7] = ori.y; nArray[8] = ori.z;
		iArray[0] = 0; iArray[1] = 1; iArray[2] = 2;
		iArray[3] = 1; iArray[4] = 0; iArray[5] = 2;
		triangles.add(new TriangleCol(ori, -Vector3f.dot(ori, P1), P1, P2, P3));
	}
	
	public RawModel addModel(Loader loader, String name, Vector3f offset, Float scale) {
		RawModel mod = ObjFileLoader.loadOBJ(name, loader);
		int[] ind = mod.getIndices();
		float[] ver = mod.getVertices();
		Vector3f P1 = new Vector3f(), P2 = new Vector3f(), P3 = new Vector3f();
		for (int i = 0; i < ind.length; i+=3) {
			P1 = new Vector3f((ver[3*ind[i]]*scale) + offset.x, (ver[3*ind[i]+1]*scale) + offset.y, (ver[3*ind[i]+2]*scale) + offset.z);
			P2 = new Vector3f((ver[3*ind[i+1]]*scale) + offset.x, (ver[3*ind[i+1]+1]*scale) + offset.y, (ver[3*ind[i+1]+2]*scale) + offset.z);
			P3 = new Vector3f((ver[3*ind[i+2]]*scale) + offset.x, (ver[3*ind[i+2]+1]*scale) + offset.y, (ver[3*ind[i+2]+2]*scale) + offset.z);
			addTriangleNoLoad(P1, P2, P3);
		}
		return ObjFileLoader.loadOBJ(name, loader);
	}
	
	public ArrayList<PlaneCol> getPlanes() {
		return planes;
	}
	
	public ArrayList<TriangleCol> getTriangles() {
		return triangles;
	}

	public ArrayList<SphereCol> getSpheres() {
		return spheres;
	}

}
