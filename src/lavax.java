import java.lang.*;
import java.util.*;
import java.io.*;
import java.nio.*;
import java.nio.file.*;

public class lavax {

  static private String  VASP_COMMAND;
  static private String  LAMMPS_COMMAND;
  static private String  INIT_POSCAR;
  static private Integer LVX_ITERATIONS;
  static private Integer VASP_NSW;
  static private double  VASP_POTIM;
  static private boolean USE_ADAPTIVE_TIMESTEP;
  static private double  MAX_DISTANCE;
  static private double  MAX_TIMESTEP;
  static private double  POTENTIAL_DEPARTURE_DISTANCE;
    
  public static class Tungsten extends Particle {

    public static double getMass() { return 183.85; }

    public Tungsten(Vec3 pos, Vec3 vel) { super(pos,vel); }
  }

  public static List<Particle> construct_bcc(double latt_const, int I, int J, int K) {
    double L = latt_const;
    List<Particle> list = new ArrayList<Particle>();

    for (int i = 0; i < I; ++i) {
      for (int j = 0; j < J; ++j) {
        for (int k = 0; k < K; ++k) {
          list.add(new Particle(new Vec3(i*L, j*L, k*L), new Vec3(0.0,0.0,0.0)));
          list.add(new Particle(new Vec3((i+0.5)*L, (j+0.5)*L, (k+0.5)*L), new Vec3(0.0,0.0,0.0)));
        }
      }
    }

    return list;
  }

  public static String make_poscar(double latt_const, Vec3 a1, Vec3 a2, Vec3 a3,
                                   List<Particle> P1, List<Particle> P2, boolean using_Cartesian) {
    StringBuilder sb = new StringBuilder();

    sb.append("BCC Xx \n");
    sb.append(latt_const + "\n");
    sb.append(a1 + "\n");
    sb.append(a2 + "\n");
    sb.append(a3 + "\n");
    sb.append(P1.size() + (P2.size() > 0 ? " " + P2.size() : "") + "\n");
    sb.append(using_Cartesian ? "Cartesian\n" : "Direct\n");

    for (Particle p : P1) {
      sb.append(p.getPos() + "\n");
    }
    for (Particle p : P2) {
      sb.append(p.getPos() + "\n");
    }
    
    sb.append("\n");
    
    for (Particle p : P1) {
      sb.append(p.getVel() + "\n");
    }
    for (Particle p : P2) {
      sb.append(p.getVel() + "\n");
    }

    return sb.toString();
  }

  public static String make_lammps_data(double xhi, double yhi, double zhi, int atom_types, 
                                        Map<Integer, Particle> P) {
    StringBuilder sb = new StringBuilder();
    
    sb.append("# W Crystal bcc #\n\n");
    sb.append(P.size() + " atoms\n");
    sb.append(atom_types + " atom types\n\n");
    sb.append(String.format("%.8f %.8f xlo xhi\n", 0.0, xhi));
    sb.append(String.format("%.8f %.8f ylo yhi\n", 0.0, yhi));
    sb.append(String.format("%.8f %.8f zlo zhi\n", 0.0, zhi));

    sb.append("\nMasses\n\n");
    sb.append(String.format("%d %.8f\n", 1, Tungsten.getMass()));
    
    sb.append("\nAtoms\n\n");

    for (Integer idx : P.keySet()) {
      // if (p instanceof Tungsten)
      sb.append(String.format("%d %d %s\n", idx, 1, P.get(idx).getPos().toString() ));
    }

    sb.append("\nVelocities\n\n");

    for (Integer idx : P.keySet()) {
      sb.append(String.format("%d %s\n", idx, P.get(idx).getVel().toString()));
    }
    
    return sb.toString();
  }
  
  public static void parse_poscar(File f, Double[] latt_const,
                                  Vec3 a1, Vec3 a2, Vec3 a3, List<Particle> P) {
    try {
      parse_poscar(new BufferedInputStream(new FileInputStream(f)), latt_const, a1, a2, a3, P);
    } catch (Throwable e) {
      System.out.println("Error " + e.getMessage());
      e.printStackTrace();
    }
  }
  
  public static void parse_poscar(InputStream is, Double[] latt_const,
                                  Vec3 a1, Vec3 a2, Vec3 a3, List<Particle> P) {
    try {
      Scanner sc = new Scanner(is);
      sc.nextLine();
      latt_const[0] = sc.nextDouble();
      
      a1.set(read_vec3(sc));
      a2.set(read_vec3(sc));
      a3.set(read_vec3(sc));

      boolean using_Direct = false;

      // Check if we are using Direct or Cartesian coordinates:
      while (sc.hasNext()) {
        String line = sc.nextLine();
        if (line.matches("C.*")) {
          using_Direct = false;
          break;
        } else if (line.matches("D.*")) {
          using_Direct = true;
          break;
        }
      }
      
      String regex = ".*[+-]?(?=[.]?[0-9])[0-9]*(?:[.][0-9]*)?(?:[Ee][+-]?[0-9]+)?.*";

      List<Vec3> pos_vec = new ArrayList<Vec3>();
      while (sc.hasNext()) {
        String line = sc.nextLine();
        if (line.matches(regex)) {
          if (using_Direct)
            pos_vec.add(transform(latt_const[0], a1, a2, a3, read_vec3(line)));
          else
            pos_vec.add((read_vec3(line)).multiply(latt_const[0]));
        } else
          break;
      }

      List<Vec3> vel_vec = new ArrayList<Vec3>();
      while (sc.hasNext()) {
        String line = sc.nextLine();
        if (line.matches(regex))
          vel_vec.add(read_vec3(line));
        else
          break;
      }

      for (int i = 0; i < pos_vec.size(); i++) {
        P.add(new Tungsten(pos_vec.get(i), vel_vec.get(i)));
      }
    }
    catch (Throwable e) {
      System.out.println("Error " + e.getMessage());
      e.printStackTrace();
    }
  }

  public static Set<Integer> parse_lammps_neighbor(File f, double cutoff) {
    Set<Integer> si = new HashSet<Integer>();
    try {
      Scanner sc = new Scanner(f);
      String regex = ".*ITEM\\: ENTRIES c_distance\\[1\\] c_neigh\\[1\\] c_neigh\\[2\\].*";
      
      while (sc.hasNext()) {
        String line = sc.nextLine();
        if (line.matches(regex)) {
          while (sc.hasNextDouble()) {
            if (sc.nextDouble() <= cutoff) {
              si.add(sc.nextInt());
              si.add(sc.nextInt());
              sc.nextLine();
            } else
              break;
          }
        }
      }
    }
    catch (Throwable e) {
      System.out.println("Error " + e.getMessage());
      e.printStackTrace();
    }
    return si;
  }

  private static Vec3 read_vec3(String line) {
    return read_vec3(new Scanner(line));
  }

  private static Vec3 read_vec3(Scanner sc) {
    double x, y, z;
    x = sc.nextDouble();
    y = sc.nextDouble();
    z = sc.nextDouble();
    
    return new Vec3(x,y,z);
  }

  public static Vec3 transform(double L, Vec3 a1, Vec3 a2, Vec3 a3, Vec3 b) {
    return new Vec3(L*(a1.getX()*b.getX() + a2.getX()*b.getX() + a3.getX()*b.getX()),
                    L*(a1.getY()*b.getY() + a2.getY()*b.getY() + a3.getY()*b.getY()),
                    L*(a1.getZ()*b.getZ() + a2.getZ()*b.getZ() + a3.getZ()*b.getZ()));
  }

  public static Vec3 PKAvel(double mass_u, Vec3 direc, double en_eV) {
    double conv_factor = 0.098237;

    double v  = conv_factor * Math.sqrt(2.0*en_eV/mass_u);
    double vc = v/Math.sqrt(direc.dot(direc));

    return direc.multiply(vc);
  }

  public static void PRINT_PROCESS(Process p) {
    PRINT_PROCESS(p.getInputStream());
  }

  public static void PRINT_PROCESS(InputStream stdin) {
    try {
      InputStreamReader isr = new InputStreamReader(stdin);
      BufferedReader br     = new BufferedReader(isr);
      String line = null;
      while ((line = br.readLine()) != null) {
        System.out.println(line);
      } 
    }
    catch (Throwable e) {
      System.out.println("Error " + e.getMessage());
      e.printStackTrace();
    }
  }

  public static void PRINT_SET(Set<Integer> s) {
    for (Integer j : s) {
      System.out.print(j + " ");
    }
    System.out.println();

  }

  public static Set<Integer> predict_neighbors() {
    try {
      // Run LAMMPS to find neighbors:
      Process p1 = Runtime.getRuntime().exec(LAMMPS_COMMAND + " -in predictor.in");
      // PRINT_PROCESS(p1);
      p1.waitFor();
      return parse_lammps_neighbor(new File("neigh.dump"), POTENTIAL_DEPARTURE_DISTANCE);
    }
    catch (Throwable e) {
      System.out.println("Error " + e.getMessage());
      e.printStackTrace();
    }
    return new HashSet<Integer>();
  }

  private static void read_conf_file(File f) {
    try {
      InputStream is = new FileInputStream(f);
      Properties prop = new Properties() {
          @Override
          public String getProperty(String key) {
            return super.getProperty(key).replaceAll("#.*", "").trim();
          }
        };

      prop.load(is);
      
      VASP_COMMAND   = prop.getProperty("VASP_COMMAND");
      LAMMPS_COMMAND = prop.getProperty("LAMMPS_COMMAND");
      INIT_POSCAR    = prop.getProperty("INIT_POSCAR");
      LVX_ITERATIONS = Integer.parseInt(prop.getProperty("LVX_ITERATIONS"));

      POTENTIAL_DEPARTURE_DISTANCE =
        Double.parseDouble(prop.getProperty("POTENTIAL_DEPARTURE_DISTANCE"));

      USE_ADAPTIVE_TIMESTEP = Boolean.parseBoolean(prop.getProperty("USE_ADAPTIVE_TIMESTEP"));
      MAX_DISTANCE   = Double.parseDouble(prop.getProperty("MAX_DISTANCE"));
      MAX_TIMESTEP   = Double.parseDouble(prop.getProperty("MAX_TIMESTEP"));
    }
    catch (Throwable e) {
      System.out.println("Error " + e.getMessage());
      e.printStackTrace();
    }
  }

  private static void init_INCAR() {
    try {
      File f = new File("INCAR");
      InputStream is = new FileInputStream(f);
      Properties prop = new Properties();

      prop.load(is);
      VASP_NSW   = Integer.parseInt(prop.getProperty("NSW"));
      VASP_POTIM = Double.parseDouble(prop.getProperty("POTIM"));
    }
    catch (Throwable e) {
      System.out.println("Error " + e.getMessage());
      e.printStackTrace();
    }
  }

  private static boolean init_check2() {
    String[] Files = new String[] {"INCAR", "KPOINTS", INIT_POSCAR, "POTCAR"};

    boolean success = true;
    for (int i = 0; i < Files.length; i++) {
      File f = new File(Files[i]);
      if (!f.exists()) {
        System.out.println("Please provide " + Files[i] + " file");
        success = false;
      }
    }
    return success;
  }

  private static void backup_files(String[] file_names, int unique_idx) {
    String folder_name = "./RUN" + unique_idx;
    Path newDirPath = Paths.get(folder_name);

    if (!Files.exists(newDirPath)) {
      try {
        Files.createDirectory(newDirPath);
        for (int i = 0; i < file_names.length; ++i) {
          Files.copy(Paths.get("./" + file_names[i]),
                     Paths.get(folder_name + "/" + file_names[i]),
                     StandardCopyOption.REPLACE_EXISTING);
        }
      }
      catch (Throwable e) {
        System.out.println("Error " + e.getMessage());
        e.printStackTrace();
      }

    }
  }

  private static void init_check1() {
    File f = new File("lavax.conf");
    if (!f.exists()) {
      System.out.println("Couldn't find file 'lavax.conf'. Creating default");
      InputStream in = lavax.class.getClassLoader().getResourceAsStream("lavax.conf");
      try {
        Files.copy(in, Paths.get("./lavax.conf"));
      }
      catch (Throwable e) {
        System.out.println("Error " + e.getMessage());
        e.printStackTrace();
      }
    }
    f = new File("predictor.in");
    if (!f.exists()) {
      System.out.println("Couldn't find file 'predictor.in'. Creating default");
      InputStream in = lavax.class.getClassLoader().getResourceAsStream("predictor.in");
      try {
        Files.copy(in, Paths.get("./predictor.in"));
      }
      catch (Throwable e) {
        System.out.println("Error " + e.getMessage());
        e.printStackTrace();
      }
    }
    System.out.println("Remember to add a LAMMPS potential file as specified in predictor.in");
    System.out.println("and a properly concatenated POTCAR file for VASP");

  }

  public static void replaceAll_in_file(File f, String regex, String replacement) {
    try {
      FileReader fr = new FileReader(f);
      BufferedReader br = new BufferedReader(fr);
      String line;

      StringBuilder sb = new StringBuilder();

      while ((line = br.readLine()) != null) {
        // s.replaceAll
        sb.append(line.replaceAll(regex, replacement) + "\n");
      }
      FileWriter fw = new FileWriter(f);
      fw.write(sb.toString());
      // System.out.println(sb.toString());
      fw.close();      
    }
    catch (Throwable e) {
      System.out.println("Error " + e.getMessage());
      e.printStackTrace();
    }
  }

  private static void update_LAMMPS_script() {
    try {
      replaceAll_in_file(new File("predictor.in"), "run.*", "run " + (VASP_NSW+3)*5);
      replaceAll_in_file(new File("predictor.in"), "^timestep.*", "timestep " + VASP_POTIM/5);
    }
    catch (Throwable e) {
      System.out.println("Error " + e.getMessage());
      e.printStackTrace();
    }
  }
  
  public static void main(String[] args) {

    init_check1();
    read_conf_file(new File("lavax.conf"));
    if (!init_check2()) {
      return;
    }
    init_INCAR();
    update_LAMMPS_script();

    // Data to be extracted from INIT_POSCAR:
    Double[] latt_const = new Double[1];
    latt_const[0] = 42.0;
    Vec3 a1 = new Vec3(0,0,0);
    Vec3 a2 = new Vec3(0,0,0);
    Vec3 a3 = new Vec3(0,0,0);
    
    List<Particle> P = new ArrayList<Particle>();
    parse_poscar(new File(INIT_POSCAR), latt_const, a1, a2, a3, P);

    // Map each particle to a unique index:
    Map<Integer, Particle> part_map = new HashMap<Integer, Particle>();

    int idx = 1;
    for (Particle p : P) {
      part_map.put(idx, p);
      ++idx;
    }
    
    // Convenience variables:
    double L = latt_const[0];
    double xrep = a1.getX();
    double yrep = a2.getY();
    double zrep = a3.getZ();

    for (int i = 0; i < LVX_ITERATIONS; i++) {
      if (USE_ADAPTIVE_TIMESTEP) {
        Particle p_max =
        Collections.max(part_map.values(),
                        new Comparator<Particle>() {
                          @Override
                          public int compare(Particle a, Particle b) {
                            Vec3 a_vel = a.getVel();
                            Vec3 b_vel = b.getVel();
                            if (a_vel.dot(a_vel) > b_vel.dot(b_vel))
                              return 1;
                            else if (a_vel.dot(a_vel) < b_vel.dot(b_vel))
                              return -1;
                            return 0;
                          }
                        });

        Vec3 max_vel = p_max.getVel();
        double speed = Math.sqrt(max_vel.dot(max_vel));
        double dt = Math.min(MAX_DISTANCE/speed, MAX_TIMESTEP);

        replaceAll_in_file(new File("INCAR"), "POTIM.*", "POTIM = " + dt);
        replaceAll_in_file(new File("predictor.in"), "^timestep.*", "timestep " + dt/5);
        replaceAll_in_file(new File("predictor.in"), "run.*", "run " + (VASP_NSW+3)*5);

        System.out.println("POTIM: " + dt);
      }
      try{
        // Write LAMMPS data file:
        PrintWriter writer = new PrintWriter("W_crystal.dat", "UTF-8");
        String content = make_lammps_data(L*xrep, L*yrep, L*zrep, 1, part_map);
        writer.print(content);
        writer.close();

        Set<Integer> si = predict_neighbors();
        System.out.println("LAMMPS run DONE");
        System.out.print("Neighbor index: ");
        PRINT_SET(si);

        // Make list of all particles that require the high precision potential:
        List<Particle> good_pot = new ArrayList<Particle>();
        for (Integer j : si) {
          good_pot.add(part_map.get(j));
        
          // Remove high precision particles from the map:
          part_map.remove(j);
        }

        // Write POSCAR:
        List<Particle> bad_pot = new ArrayList<Particle>(part_map.values());
        writer = new PrintWriter("POSCAR", "UTF-8");
        content = make_poscar(1.0,
                              a1.multiply(latt_const[0]),
                              a2.multiply(latt_const[0]),
                              a3.multiply(latt_const[0]),
                              bad_pot, good_pot, true);
        writer.print(content);
        writer.close();
        
        // Run VASP:
        Process p1 = Runtime.getRuntime().exec(VASP_COMMAND);
        PRINT_PROCESS(p1.getInputStream());
        PRINT_PROCESS(p1.getErrorStream());
        p1.waitFor();
        System.out.println("VASP run DONE");
        
        // Parse CONTCAR:
        P.clear();
        parse_poscar(new File("CONTCAR"), latt_const, a1, a2, a3, P);

        // Construct new Particle map:
        part_map.clear();
        idx = 1;
        for (Particle p : P) {
          part_map.put(idx, p);
          ++idx;
        }
        
        String[] files = {"XDATCAR", "CONTCAR", "CHG",
                          "CHGCAR", "DOSCAR", "EIGENVAL",
                          "OSZICAR", "PCDAT", "vasprun.xml",
                          "OUTCAR", "INCAR", "WAVECAR",
                          "IBZKPT", "POSCAR"};
        backup_files(files, i);
        System.out.println("------------------------------");
        
      } catch (Throwable e) {
        System.out.println("Error " + e.getMessage());
        e.printStackTrace();
      }
    }    
  }
}
