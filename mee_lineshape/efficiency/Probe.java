import org.jlab.jnp.hipo4.io.HipoReader;
import org.jlab.jnp.hipo4.data.*;
public class Probe {
  public static void main(String[] a) {
    HipoReader r = new HipoReader(); r.open(a[0]);
    String[] names = {"MC::Particle","MC::Lund","MC::RecMatch","MC::GenMatch","REC::Particle","REC::Calorimeter","RUN::config"};
    Bank[] b = new Bank[names.length]; int[] c = new int[names.length];
    for (int i = 0; i < names.length; i++) if (r.getSchemaFactory().hasSchema(names[i])) b[i] = r.getBank(names[i]);
    Event ev = new Event(); int n = 0;
    while (r.hasNext() && n < 200) { r.nextEvent(ev); n++;
      for (int i = 0; i < names.length; i++) if (b[i] != null) { ev.read(b[i]); if (b[i].getRows() > 0) c[i]++; } }
    System.out.println("events read " + n + " of " + r.getEventCount());
    for (int i = 0; i < names.length; i++) System.out.println("  " + names[i] + (b[i]==null ? "  (no schema)" : "  non-empty in " + c[i]));
  }
}
