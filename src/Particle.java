import java.util.List;

public class Particle {

  Vec3 pos;
  Vec3 vel;
  List<Integer> neighbors;
  
  public Particle(Vec3 pos, Vec3 vel) {
    this.pos = pos;
    this.vel = vel;
    this.neighbors = null;
  }
  
  public Particle(Vec3 pos, Vec3 vel, List<Integer> neighbors) {
    this(pos, vel);
    this.neighbors = neighbors;
  }

  public Vec3 getVel() {
    return vel;
  }
  
  public Vec3 getPos() {
    return pos;
  }

  public List<Integer> getNeighbors() {
    return neighbors;
  }

  public void update(Vec3 pos, Vec3 vel, List<Integer> neighbors) {
    this.pos = pos;
    this.vel = vel;
    this.neighbors = neighbors;
  }

  public void update(Vec3 pos, Vec3 vel) {
    this.pos = pos;
    this.vel = vel;
  }

  @Override
  public String toString() {
    return pos + "\n" + vel;
  }

  // public static class Tungsten extends Particle {

  //   public Tungsten(Vec3 pos, Vec3 vel, List<Integer> neighbors) {
  //     super(pos, vel, neighbors);
  //   }

  //   public static double getMass() {
  //     return 183.85;
  //   }
  // }
}
