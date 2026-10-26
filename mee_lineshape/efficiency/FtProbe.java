import org.jlab.jnp.hipo4.io.HipoReader;
import org.jlab.jnp.hipo4.data.*;
// events with generated scattered e' (first pid 11 in MC::Particle) at p 0.5-0.8 GeV, theta 3-5 deg:
// how many have an FTCAL entry in REC::ForwardTagger, its energy, and an REC::Particle in FT.
public class FtProbe { public static void main(String[] a) {
  HipoReader r = new HipoReader(); r.open(a[0]);
  String[] want = {"REC::ForwardTagger","FTCAL::clusters","FT::particles","FTCAL::hits","FTHODO::hits","REC::Particle","MC::Particle"};
  for (String w : want) System.out.println("  schema " + w + ": " + r.getSchemaFactory().hasSchema(w));
  Bank mc = r.getBank("MC::Particle"), ft = r.getBank("REC::ForwardTagger"), rp = r.getBank("REC::Particle");
  Bank cl = r.getSchemaFactory().hasSchema("FTCAL::clusters") ? r.getBank("FTCAL::clusters") : null;
  Event ev = new Event(); int nsel = 0, nft = 0, nrecft = 0, ncl = 0; double esum = 0; int ne = 0;
  double[] eb = new double[12];
  while (r.hasNext()) { r.nextEvent(ev); ev.read(mc); ev.read(ft); ev.read(rp); if (cl != null) ev.read(cl);
    double p = -1, th = -1;
    for (int i = 0; i < mc.getRows(); i++) if (mc.getInt("pid", i) == 11) { double px=mc.getFloat("px",i),py=mc.getFloat("py",i),pz=mc.getFloat("pz",i);
      p = Math.sqrt(px*px+py*py+pz*pz); th = Math.toDegrees(Math.atan2(Math.hypot(px,py),pz)); break; }
    if (!(p > 0.5 && p < 0.8 && th > 3 && th < 5)) continue;
    nsel++;
    boolean hasft = false;
    for (int i = 0; i < ft.getRows(); i++) if (ft.getByte("detector", i) == 10) { hasft = true; double e = ft.getFloat("energy", i); esum += e; ne++; int k = (int)Math.min(11, e / 0.1); eb[k]++; }
    if (hasft) nft++;
    if (cl != null && cl.getRows() > 0) ncl++;
    for (int i = 0; i < rp.getRows(); i++) { int s = Math.abs(rp.getShort("status", i)); if (s >= 1000 && s < 2000 && rp.getInt("pid", i) == 11) { nrecft++; break; } }
  }
  System.out.printf("selected %d  with FTCAL in REC::ForwardTagger %d (%.3f)  FTCAL::clusters non-empty %d  REC e- in FT %d (%.3f)  <E_FTCAL> %.3f%n",
     nsel, nft, nft/(double)nsel, ncl, nrecft, nrecft/(double)nsel, ne > 0 ? esum/ne : -1);
  StringBuilder sb = new StringBuilder("  E_FTCAL histogram (0.1 GeV bins): "); for (double x : eb) sb.append((int)x).append(' '); System.out.println(sb);
} }
