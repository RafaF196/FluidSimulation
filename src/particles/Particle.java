package particles;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import java.util.Random;
import java.util.Vector;

import org.lwjgl.util.vector.Vector3f;

import collisions.CollisionMaster;
import collisions.PlaneCol;
import collisions.SphereCol;
import collisions.TriangleCol;

public class Particle {
	
	private Vector3f position;
	private Vector3f offset = new Vector3f(0, 0, 0);
	private Vector3f oldoffset = new Vector3f(0, 0, 0);
	private Vector3f velocity;
	private Vector3f oldvelocity = new Vector3f(0, 0, 0);
	private float mass;
	private float gravity;
	private float lifeLength;
	
	private float reference_density;
	private float density;
	private float pressure;
	
	private static float c_s = 5f;
	private static float dyn_vis = 0.001f;
	
	Vector<Integer> neigh_indexes = new Vector<Integer>();
	Vector<Float> neigh_dist = new Vector<Float>();
	
	private float elapsedTime = 0;
	
	private Boolean isFixed = false;
	
	Vector3f auxVec3 = new Vector3f();
	Vector3f auxVec3b = new Vector3f();
	Vector3f auxVec3c = new Vector3f();
	
	Vector3f newOffset = new Vector3f();
	Vector3f newVelocity = new Vector3f();
	Vector3f diff = new Vector3f();
	Vector3f force = new Vector3f();
	
	Vector3f v = new Vector3f();
	Vector3f pa = new Vector3f();
	Vector3f norm = new Vector3f();
	Vector3f normVel = new Vector3f();
	Vector3f tanVel = new Vector3f();
	
	Vector3f correctPos = new Vector3f();
	Vector3f correctVel = new Vector3f();
	
	Vector3f P1 = new Vector3f();
	Vector3f P2 = new Vector3f();
	Vector3f v1 = new Vector3f();
	Vector3f v2 = new Vector3f();
	
	public Particle(Vector3f position, Vector3f offset, Integer numParticles, Vector3f velocity, float mass, float reference_density,
			float delta, float gravity, float lifeLength) {
		this.position = position;
		this.velocity = velocity;
		this.mass = mass;
		this.reference_density = reference_density;
		
		this.offset.set(offset);
		// Previous position
		this.oldoffset.set(velocity);
		this.oldoffset.scale(-delta);
		Vector3f.add(oldoffset, offset, oldoffset);

		this.gravity = gravity;
		this.lifeLength = lifeLength;
	}

	private void setOffset(Vector3f pos) {
		Vector3f.sub(pos, this.position, this.offset);
	}
	
	private void setOldOffset(Vector3f pos) {
		Vector3f.sub(pos, this.position, this.oldoffset);
	}
	
	protected Vector3f getPosition() {
		return Vector3f.add(position, offset, null);
	}
	
	private Vector3f getPrevPosition() {
		return Vector3f.add(position, oldoffset, null);
	}
	
	public Vector3f getOffset() {
		return offset;
	}

	public Vector3f getVelocity() {
		return velocity;
	}
	
	public Vector3f getOldVelocity() {
		return oldvelocity;
	}
	
	public Float getDensity() {
		return density;
	}
	
	public Float getPressure() {
		return pressure;
	}
	
	private Vector3f getSegmentPlaneIntersection(Vector3f n, Float d) {
		Float v1 = Vector3f.dot(n, getPrevPosition());
		Float v2 = Vector3f.dot(n, getPosition());
		Vector3f.sub(offset, oldoffset, auxVec3);
		Vector3f P = new Vector3f();
		auxVec3.scale(-(-d-v1)/(v1-v2));
		Vector3f.add(getPrevPosition(), auxVec3, P);
		return P;
	}
	
	private Float calculateArenaTriangle(Vector3f Vi, Vector3f Vj, Vector3f Vk) {
		Vector3f.sub(Vj, Vi, auxVec3);
		Vector3f.sub(Vk, Vi, auxVec3b);
		Vector3f.cross(auxVec3, auxVec3b, auxVec3c);
		return auxVec3c.length()/2;
	}
	
	private Boolean isIntersectingTriangle(Vector3f V1, Vector3f V2, Vector3f V3, Vector3f X) {
		return calculateArenaTriangle(X,V2,V3) + calculateArenaTriangle(V1,X,V3) + calculateArenaTriangle(V1,V2,X) -
				calculateArenaTriangle(V1,V2,V3) < 0.001f; // Small offset to avoid small errors
	}
	
	private Vector3f getSegmentSphereIntersection(Vector3f C, Float R) {
		Vector3f v1 = new Vector3f(C.x - getPrevPosition().x, C.y - getPrevPosition().y, C.z - getPrevPosition().z);
		Vector3f v2 = new Vector3f(C.x - getPosition().x, C.y - getPosition().y, C.z - getPosition().z);
		Float r1 = v1.length();
		Float r2 = v2.length();
		if ((r1 > R && r2 < R) || (r1 < R && r2 > R) || (r2 < R && r1 < R)) {
			Vector3f.sub(offset, oldoffset, auxVec3);
			Vector3f P = new Vector3f();
			if (r1-r2 == 0) auxVec3.scale(0);
			else auxVec3.scale(-(R-r1)/(r1-r2));
			Vector3f.add(getPrevPosition(), auxVec3, P);
			return P;
		}
		return null;
	}

	private void updateNeighbors(ArrayList<Particle> particles, ArrayList<ArrayList<ArrayList<Integer>>> particle_grid,
			Integer index) {
		this.neigh_indexes = new Vector<Integer>();
		this.neigh_dist = new Vector<Float>();
		Vector3f auxV4, auxV5 = this.getPosition(), difference = new Vector3f(0,0,0);
		ArrayList<Integer> auxArr;
		int xValue = (int) ((particles.get(index).getPosition().x + 1.1f) / 0.2f);
		int zValue = (int) ((particles.get(index).getPosition().z + 1.1f) / 0.2f);
		for (int i = xValue-1; i <= xValue+1; i++) {
			for (int j = zValue-1; j <= zValue+1; j++) {
				auxArr = particle_grid.get(i).get(j);
				for (int k = 0; k < auxArr.size(); k++) {
					auxV4 = particles.get(auxArr.get(k)).getPosition();
					Vector3f.sub(auxV4, auxV5, difference);
					if (difference.length() < 0.2f) {
						neigh_indexes.add(auxArr.get(k));
						neigh_dist.add(difference.length());
					}
				}
			}
		}
	}

	private void computeDensity(ArrayList<Particle> particles, Integer index) {
		float temp_dens = 0f;
		for (int i = 0; i < neigh_indexes.size(); i++) {
			temp_dens += this.mass * (0.2f - neigh_dist.get(i))/0.2;
		}
		this.density = temp_dens;
	}
	
	private void computePressure(ArrayList<Particle> particles, Integer index) {
		this.pressure = c_s * c_s * (this.density - this.reference_density);
	}
	
	private Vector3f getPressureForce(ArrayList<Particle> particles, Integer index) {
		Vector3f temp_PF = new Vector3f(0,0,0), r_ij = new Vector3f(0,0,0);
		Float P_ij = 0.0f, d_j, p_j;
		for (int i = 0; i < this.neigh_indexes.size(); i++) {
			d_j = particles.get(neigh_indexes.get(i)).getDensity();
			p_j = particles.get(neigh_indexes.get(i)).getPressure();
			P_ij = -mass*((pressure/(density*density)) + (p_j/(d_j*d_j)));
			Vector3f.sub(getPosition(), particles.get(neigh_indexes.get(i)).getPosition(), r_ij);
			P_ij *= (float) ((0.2f - neigh_dist.get(i))/200f);
			Vector3f.add(temp_PF, new Vector3f(mass * P_ij * r_ij.x, mass * P_ij * r_ij.y, mass * P_ij * r_ij.z), temp_PF);
		}
		
		return temp_PF;
	}

	private Vector3f getViscosityForce(ArrayList<Particle> particles, Integer index) {
		Vector3f temp_VF = new Vector3f(0,0,0), r_ij = new Vector3f(0,0,0);
		Float V_ij = 0.0f, d_j;
		for (int i = 0; i < this.neigh_indexes.size(); i++) {
			d_j = particles.get(neigh_indexes.get(i)).getDensity();
			V_ij = dyn_vis*mass*(((dyn_vis/density) - (dyn_vis/d_j))/(density*d_j));
			Vector3f.sub(getPosition(), particles.get(neigh_indexes.get(i)).getPosition(), r_ij);
			V_ij *= (float) ((0.2f - neigh_dist.get(i))/50f);
			Vector3f.add(temp_VF, new Vector3f(mass * V_ij * r_ij.x, mass * V_ij * r_ij.y, mass * V_ij * r_ij.z), temp_VF);
		}
		
		return temp_VF;
	}
	
	public void updateFluidInfo(ArrayList<Particle> particles, ArrayList<ArrayList<ArrayList<Integer>>> particle_grid,
			Integer index) {
		// Fluid calculations
		updateNeighbors(particles, particle_grid, index);
		computeDensity(particles, index);
		computePressure(particles, index);
	}
	
	public void update(ArrayList<Particle> particles, Float L, Integer index, Float delta, CollisionMaster cm,
			Float friction, Float bouncing, Boolean interactive_force, Integer numParticles){
		
		Vector3f gravityForce = new Vector3f(0, -1*mass*gravity, 0);
		force.set(gravityForce);
		// Add a_pressure, a_viscosity
		Vector3f.add(force, getPressureForce(particles, index), force);
		Vector3f.add(force, getViscosityForce(particles, index), force);
		if (interactive_force) Vector3f.add(force, new Vector3f(0, -1*mass*20, 0), force);
		force.scale(delta*delta/mass);
	}
	
	public boolean update2(ArrayList<Particle> particles, Float L, Integer index, Float delta, CollisionMaster cm,
			Float friction, Float bouncing, Integer numParticles){
		elapsedTime += delta;
		if (this.isFixed) return elapsedTime > lifeLength;
		
		// Position
		newOffset.set(offset);
		Vector3f.sub(offset, oldoffset, diff);
		diff.scale(0.999f);
		Vector3f.add(newOffset, diff, newOffset);
		
		Vector3f.add(newOffset, force, newOffset);
		
		oldoffset.set(offset);
		offset.set(newOffset);
		
		// Velocity
		oldvelocity.set(velocity);
		Vector3f.sub(offset, oldoffset, auxVec3);
		velocity = new Vector3f(auxVec3.x / delta, auxVec3.y / delta, auxVec3.z / delta);
		
		checkCollision(delta, cm, friction, bouncing);

		return elapsedTime > lifeLength;
	}
	
	private void handleCollision(Vector3f normal, Float d, Float delta, Float friction, Float bouncing) {
		
		// Position
		Vector3f p_prime = new Vector3f();
		p_prime.set(this.getPosition());
		Float aux = (1+bouncing)*(Vector3f.dot(normal, p_prime) + d);
		auxVec3.set(normal);
		auxVec3.scale(aux);
		Vector3f.sub(p_prime, auxVec3, correctPos);
		
		// Velocity
		Vector3f v_prime = new Vector3f();
		v_prime.set(this.velocity);
		aux = (1+bouncing)*(Vector3f.dot(normal, v_prime));
		auxVec3.set(normal);
		auxVec3.scale(aux);
		Vector3f.sub(v_prime, auxVec3, correctVel);
		
		float fvn = Vector3f.dot(normal, correctVel);
		normVel.set(normal);
		normVel.scale(fvn);
		Vector3f.sub(correctVel, normVel, tanVel);
		tanVel.scale(friction);
		Vector3f.sub(correctVel, tanVel, correctVel);
		
		// Previous position
		auxVec3.set(correctVel);
		auxVec3.scale(-delta);
		Vector3f.add(correctPos, auxVec3, auxVec3b);
		
		this.setOffset(correctPos);
		this.setOldOffset(auxVec3b);
		this.velocity.set(correctVel);
	}
	
	private void checkCollision(Float delta, CollisionMaster cm, Float friction, Float bouncing) {
		Vector3f ptdt_prime = getPosition();
		Vector3f pt = getPrevPosition();
		
		// Triangle collisions
		ArrayList<TriangleCol> triangles = cm.getTriangles();		
		for (int i = 0; i < triangles.size(); i++) {
			TriangleCol t = triangles.get(i);
			norm.set(t.getNormal());
			Float d = t.getD();
			Float val = (Vector3f.dot(norm, ptdt_prime) + d)*(Vector3f.dot(norm, pt) + d);
			if (val <= -0) {
				Vector3f P = getSegmentPlaneIntersection(norm, d);
				if (isIntersectingTriangle(t.getP1(), t.getP2(), t.getP3(), P)) handleCollision(norm, t.getD(), delta, friction, bouncing);
			}
		}
		
		// Sphere collisions
		ArrayList<SphereCol> spheres = cm.getSpheres();
		for (int i = 0; i < spheres.size(); i++) {
			SphereCol s = spheres.get(i);
			Vector3f C = s.getCenter();
			Vector3f P = getSegmentSphereIntersection(C, s.getRadius());
			if (P != null) {
				Vector3f normal = new Vector3f(P.x - C.x, P.y - C.y, P.z - C.z);
				normal.normalise();
				handleCollision(normal, -Vector3f.dot(normal, P), delta, friction, bouncing);
			}
		}
		
		// Plane collisions
		ArrayList<PlaneCol> planes = cm.getPlanes();
		for (int i = 0; i < planes.size(); i++) {
			PlaneCol p = planes.get(i);
			norm.set(p.getNormal());
			Float d = p.getValue();
			//Float val = (Vector3f.dot(norm, ptdt_prime) + d)*(Vector3f.dot(norm, pt) + d);
			//if (val <= 0/* || this.getPosition().y < -d*/) handleCollision(norm, p.getValue(), delta, friction, bouncing);
			if (norm.y > 0.1f && this.getPosition().y < -d) handleCollision(norm, p.getValue(), delta, friction, bouncing);
			if (norm.x > 0.1f && d > 0 && this.getPosition().x < -d) handleCollision(norm, p.getValue(), delta, friction, bouncing);
			if (norm.x > 0.1f && d < 0 && this.getPosition().x > -d) handleCollision(norm, p.getValue(), delta, friction, bouncing);
			if (norm.z > 0.1f && d > 0 && this.getPosition().z < -d) handleCollision(norm, p.getValue(), delta, friction, bouncing);
			if (norm.z > 0.1f && d < 0 && this.getPosition().z > -d) handleCollision(norm, p.getValue(), delta, friction, bouncing);
			
		}
	}

	public void unfix() {
		this.isFixed = false;
	}

}
