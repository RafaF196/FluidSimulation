package models;

public class RawModel {
	
	private String modelName;
	
	private int vaoID; // ID of the VAO
	private int vertexCount; // number of vertex of the model
	
	private float[] verticesArray; // array of all the vertex positions
	private float[] texturesArray; // array of all the texture coordinates
	private float[] normalsArray; // array of all the normal vectors of the vertices
	private int[] indicesArray; // indices of vertices to obtain the faces
	
	public RawModel(String modelName, int vaoID, int vertexCount, float[] verticesArray, float[] texturesArray,
			float[] normalsArray, int[] indicesArray) {
		this.modelName = modelName;
		this.vaoID = vaoID;
		this.vertexCount = vertexCount;
		
		this.verticesArray = verticesArray;
		this.texturesArray = texturesArray;
		this.normalsArray = normalsArray;
		this.indicesArray = indicesArray;
		
	}

	public String getModelName() {
		return modelName;
	}
	
	public int getVaoID() {
		return vaoID;
	}

	public int getVertexCount() {
		return vertexCount;
	}
	
	public float[] getVertices() {
		return verticesArray;
	}
	
	public float[] getTextures() {
		return texturesArray;
	}
	
	public float[] getNormals() {
		return normalsArray;
	}
	
	public int[] getIndices() {
		return indicesArray;
	}
	
}
