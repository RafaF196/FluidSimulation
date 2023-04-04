package particles;

import java.util.ArrayList;
import java.util.Random;
import java.util.Set;

import org.lwjgl.util.vector.Vector3f;

import collisions.CollisionMaster;
import models.TexturedModel;

public class ParticleSystem {
	
	private TexturedModel model;
	private Vector3f initialPos;
	private Integer numParticles;
	private Float L;
	
	ArrayList<Particle> particles = new ArrayList<Particle>();
	
	// List of particle indexes in any cell of the grid (12x12 including borders, 0.3 cell size)
	ArrayList<ArrayList<ArrayList<Integer>>> particle_grid = new ArrayList<ArrayList<ArrayList<Integer>>>();
	
	Random rand = new Random();
	
	public ParticleSystem(TexturedModel model, Vector3f initialPos, Float separation, Integer numParticles) {
		this.model = model;
		this.initialPos = initialPos;
		this.L = separation;
		this.numParticles = numParticles;
		addFluid(initialPos, separation, numParticles);
	}

	public TexturedModel getModel() {
		return model;
	}
	
	public Vector3f getPosition() {
		return initialPos;
	}
	
	public Integer getNumParticles() {
		return numParticles;
	}
	
	public Float getL() {
		return L;
	}
	
	public ArrayList<Vector3f> getOffsets() {
		ArrayList<Vector3f> os = new ArrayList<Vector3f>();
		for (int i = 0; i < particles.size(); i++) os.add(particles.get(i).getOffset());
		return os;
	}
	
	public void addFluid(Vector3f initPos, Float separation, Integer numParticles) {
		Vector3f totalOffset = new Vector3f(0,0,0);
		Float reference_density = 1000f; // may change
		Random rand = new Random();
		for (int i = 0; i < numParticles; i++) {
			for (int j = 0; j < numParticles; j++) {
				for (int k = 0; k < numParticles; k++) {
					Vector3f initV = new Vector3f(0.0f, 0.0f, 0.0f);
					totalOffset = new Vector3f(separation*i/* + rand.nextFloat()/100*/,
							separation*k,
							separation*j);
					Particle p = new Particle(initialPos, totalOffset, numParticles, initV, 0.01f,
							reference_density, 0.01f, 9.8f, 1000);
					particles.add(p);
				}
			}
		}
	} 
	
	private void fillGrid() {
		particle_grid = new ArrayList<>(11);
		for (int i = 0; i < 11; i++) {
			particle_grid.add(new ArrayList<ArrayList<Integer>>(11));
		    for (int j = 0; j < 11; j++) {
		    	particle_grid.get(i).add(new ArrayList<Integer>(11));
		    }
		}
		
		for (int i = 0; i < numParticles*numParticles*numParticles; i++) {
			int xValue = (int) ((particles.get(i).getPosition().x + 1.1f) / 0.2f);
			int zValue = (int) ((particles.get(i).getPosition().z + 1.1f) / 0.2f);
			this.particle_grid.get(xValue).get(zValue).add(i);
		}
		
		/*
		System.out.println("==================================================================");
		for (int i = 0; i < 11; i++) {
			System.out.println(this.particle_grid.get(i).get(0).size() + " " + 
					this.particle_grid.get(i).get(1).size() + " " + 
					this.particle_grid.get(i).get(2).size() + " " + 
					this.particle_grid.get(i).get(3).size() + " " + 
					this.particle_grid.get(i).get(4).size() + " " + 
					this.particle_grid.get(i).get(5).size() + " " + 
					this.particle_grid.get(i).get(6).size() + " " + 
					this.particle_grid.get(i).get(7).size() + " " + 
					this.particle_grid.get(i).get(8).size() + " " + 
					this.particle_grid.get(i).get(9).size() + " " + 
					this.particle_grid.get(i).get(10).size());
		}
		System.out.println("==================================================================");
		*/
	}
	
	public void update(Float delta, CollisionMaster cm, Float friction, Float bouncing, Boolean interactive_force) {
		Boolean delete;
		fillGrid();
		for (int i = 0; i < particles.size(); i++) { // Calculates density and pressure
			particles.get(i).updateFluidInfo(particles, particle_grid, i);
		}
		for (int i = 0; i < particles.size(); i++) { // Calculates the new forces
			particles.get(i).update(particles, this.L, i, delta, cm, friction, bouncing, interactive_force, numParticles);
		}
		for (int i = 0; i < particles.size(); i++) { // Applies the forces
			delete = particles.get(i).update2(particles, this.L, i, delta, cm, friction, bouncing, numParticles);
			if (delete) {
				particles.remove(i);
				i--;
			}
		}
	}
	
	public void unfixParticles() {
		for (int i = 0; i < particles.size(); i++) {
			particles.get(i).unfix();
		}
	}
	
}
