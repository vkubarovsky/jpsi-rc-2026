import org.jlab.jnp.hipo4.io.HipoReader;
public class Ver { public static void main(String[] a) { HipoReader r = new HipoReader(); r.open(a[0]);
  try { for (java.lang.reflect.Method m : r.getClass().getMethods()) if (m.getName().toLowerCase().contains("config") && m.getParameterCount()==0) System.out.println(m.getName()+": "+m.invoke(r)); } catch (Exception e) { System.out.println(e); } } }
