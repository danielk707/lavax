
public class Vec3 {
  double x;
  double y;
  double z;

  public double getX() {
    return x;
  }

  public double getY() {
    return y;
  }

  public double getZ() {
    return z;
  }

  public Vec3(double x, double y, double z) {
    this.x = x;
    this.y = y;
    this.z = z;
  }

  public void set(Vec3 v) {
    this.x = v.getX();
    this.y = v.getY();
    this.z = v.getZ();
  }

  public Vec3 add(Vec3 v) {
    return new Vec3(this.x + v.getX(), this.y + v.getY(), this.z + v.getZ());
  }

  public Vec3 subtract(Vec3 v){
    return new Vec3(this.x- v.getX(), this.y - v.getY(), this.z - v.getZ());
  }

  public Vec3 multiply(double a){
    return new Vec3(a*this.x,a*this.y,a*this.z);
  }

  public double dot(Vec3 v){
    return v.getX()*this.x+this.y*v.getY()+this.z*v.getZ();
  }

  @Override
  public String toString() {
    return String.format("%.8f %.8f %.8f", x, y, z);
  }

  @Override
  public boolean equals(Object obj) {
    return toString().equals(obj.toString());
  }
}
