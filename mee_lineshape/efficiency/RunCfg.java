import org.jlab.jnp.hipo4.io.HipoReader;
import org.jlab.jnp.hipo4.data.*;
public class RunCfg { public static void main(String[] a) { HipoReader r = new HipoReader(); r.open(a[0]);
  Bank b = r.getBank("RUN::config"); Event ev = new Event(); r.nextEvent(ev); ev.read(b);
  System.out.println("run " + b.getInt("run",0) + "  torus " + b.getFloat("torus",0) + "  solenoid " + b.getFloat("solenoid",0)); } }
