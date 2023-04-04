package particles;

import java.util.ArrayList;

import org.lwjgl.util.vector.Matrix4f;
import org.lwjgl.util.vector.Vector3f;

import toolbox.MatrixOps;

import entities.Camera;
import shaders.ShaderProgram;

public class ParticleShader extends ShaderProgram {
	
	private static final String VERTEX_FILE = "/particles/particleVertexShader.txt";
	private static final String FRAGMENT_FILE = "/particles/particleFragmentShader.txt";
	
	private int location_transformationMatrix;
	private int location_projectionMatrix;
	private int location_viewMatrix;
	private int location_modelTexture;
	private int location_cameraPosition;
	private int location_offsetPos[];

	public ParticleShader() {
		super(VERTEX_FILE, FRAGMENT_FILE);
	}

	@Override
	protected void bindAttributes() {
		super.bindAttribute(0, "position");
		super.bindAttribute(1, "textureCoordinates");
		super.bindAttribute(2, "normal");
	}

	@Override
	protected void getAllUniformLocations() {
		location_transformationMatrix = super.getUniformLocation("transformationMatrix");
		location_projectionMatrix = super.getUniformLocation("projectionMatrix");
		location_viewMatrix = super.getUniformLocation("viewMatrix");
		location_modelTexture = super.getUniformLocation("modelTexture");
		location_cameraPosition = super.getUniformLocation("cameraPosition");
		
		location_offsetPos = new int[2000];
		for(int i = 0; i < 2000; i++){
			location_offsetPos[i] = super.getUniformLocation("offsetPos[" + i + "]");
		}
	}
	
	protected void connectTextureUnits() {
		super.loadInt(location_modelTexture, 0);
	}
	
	public void loadTransformationMatrix(Matrix4f matrix) {
		super.loadMatrix(location_transformationMatrix, matrix);
	}
	
	public void loadViewMatrix(Camera camera) {
		Matrix4f viewMatrix = MatrixOps.createViewMatrix(camera);
		super.loadMatrix(location_viewMatrix, viewMatrix);
		super.loadVector(location_cameraPosition, camera.getPosition());
	}
	
	public void loadProjectionMatrix(Matrix4f projection) {
		super.loadMatrix(location_projectionMatrix, projection);
	}
	
	// Generates the offset for each model given the amount of elements in one row of the grid
	public void loadOffsetPosition(ArrayList<Vector3f> os) {
		for (int i = 0; i < os.size(); i++) {
			super.loadVector(location_offsetPos[i], os.get(i));
		}
	}

}
