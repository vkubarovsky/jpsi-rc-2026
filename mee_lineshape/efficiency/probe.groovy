import org.jlab.jnp.hipo4.io.HipoReader
import org.jlab.jnp.hipo4.data.Event
def r = new HipoReader(); r.open(args[0])
def ev = new Event()
def names = ["MC::Particle","MC::Lund","MC::RecMatch","MC::GenMatch","REC::Particle","RUN::config"]
def banks = names.collectEntries { n -> [(n): (r.getSchemaFactory().hasSchema(n) ? r.getBank(n) : null)] }
int n = 0; def cnt = [:].withDefault{0}
while (r.hasNext() && n < 200) { r.nextEvent(ev); n++
  banks.each { k, b -> if (b != null) { ev.read(b); if (b.getRows() > 0) cnt[k]++ } } }
println "events read: $n  (file has ${r.getEventCount()})"
cnt.each { k, v -> println "  $k non-empty in $v events" }
