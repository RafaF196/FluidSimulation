package entityRenderer;

import org.lwjgl.util.vector.Matrix4f;
import org.lwjgl.util.vector.Vector2f;

import toolbox.MatrixOps;

import entities.Camera;
import shaders.ShaderProgram;

public class EntityShader extends ShaderProgram {
	
	private static final String VERTEX_FILE = "/entityRenderer/entityVertexShader.txt";
	private static final String FRAGMENT_FILE = "/entityRenderer/entityFragmentShader.txt";
	
	private int location_transformationMatrix;
	private int location_projectionMatrix;
	private int location_viewMatrix;
	private int location_modelTexture;
	private int location_cameraPosition;
	private int location_offsetPos[];

	public EntityShader() {
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
		
		location_offsetPos = new int[225];
		for(int i = 0; i < 225; i++){
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
	public void loadOffsetPosition(Integer rows) {
		Integer index = 0;
		for (int i = -rows; i < rows; i+=2) {
			for (int j = -rows; j < rows; j+=2) {
				Vector2f pos = new Vector2f(i*2 +2, j*2 +2);
				super.load2DVector(location_offsetPos[index], pos);
				index++;
			}
		}
	}

}
