<script>
export default {
	data: function() {
		return {
			errormsg: null,
			loading: false,
			some_data: null,
		}
	},
	methods: {
        data(){
        },
	},
}
</script>

<template>
	<div>
        <h1 id="microalgae-parsing-searching">MicroAlgae Parsing &amp;
Searching</h1>
<p>Credits <a
href="https://github.com/federico-rosatelli"><code>federico-rosatelli</code></a>
<a href="https://github.com/AxnNxs"><code> Mat</code></a> <a
href="https://github.com/Loriv3"><code>Loriv3</code></a> <a
href="https://github.com/BboyCaligola"><code>Calli</code></a> <a
href="https://github.com/Samseys"><code>Samsey</code></a></p>
<h1 id="requisiti">Requisiti:</h1>
<p>Anaconda: “https://www.anaconda.com”</p>
<p>Biopython: “https://biopython.org/wiki/Documentation”</p>
<p>Hashlib: “https://docs.python.org/3/library/hashlib.html”</p>
<p>Argparse: “https://docs.python.org/3/library/argparse.html”</p>
<p>Pymongo: “https://pymongo.readthedocs.io/en/stable/”</p>
<p>MongoDB CE “Jammy”:
“https://www.mongodb.com/docs/manual/tutorial/install-mongodb-on-ubuntu/”</p>
<p>NVM:
“https://www.linode.com/docs/guides/how-to-install-use-node-version-manager-nvm/”</p>
<p>Nodejs: “https://github.com/nodejs/help/wiki/Installation”</p>
<p>NPM:
“https://docs.npmjs.com/downloading-and-installing-node-js-and-npm”</p>
<p>Nodejs Express: “http://expressjs.com/”</p>
<p>GoLang: “https://go.dev/doc/install”</p>
<p>VUEjs: “https://cli.vuejs.org/guide/installation.html”</p>
<p>Docker: “https://docs.docker.com/engine/install/ubuntu/”</p>
<p>Docker compose:
“https://docs.docker.com/desktop/install/linux-install/” (Integrato in
Desktop)</p>
<h1 id="fonti-ufficiali">Fonti ufficiali:</h1>
<p>Shigen: “https://shigen.nig.ac.jp” DB giapponese</p>
<p>NCBI: “https://www.ncbi.nlm.nih.gov/” DB americano</p>
<p>La maggior parte dei dati presenti localmente nel Database interno
sono stati scaricati manualmente a causa delle restrizioni sui
download.</p>
<h1 id="altri-link">Altri link:</h1>
<p>Drive:
“https://drive.google.com/drive/folders/19RXRHEb-7-O9gaUjXz5ho-Q2_HsbKlEW”</p>
<p>Github: “https://github.com/federico-rosatelli/biologia”</p>
<h1 id="guida-al-primo-utilizzo">Guida al primo utilizzo</h1>
<p>Per la prima installazione su un qualsiasi PC, seguire i seguenti
passaggi (si raccomanda Ubuntu):</p>
<ul>
<li><p>Scaricare il progetto tramite git, o copiare il progetto
all’indirizzo inserito tra “Altri link” sopra questo paragrafo.</p></li>
<li><p>Installare tutti i requisiti, settando le variabili di ambiente a
livello globale</p></li>
<li><p>Inizializzare MongoDB: Dare i seguenti comandi da terminale:</p>
<p><code>sudo chown -R mongodb:mongodb /var/lib/mongodb</code></p>
<p><code>sudo chown mongodb:mongodb /tmp/mongodb-27017.sock</code></p>
<p>dopodiché:</p>
<p><code>sudo systemctl start mongod</code></p>
<p>per rendere automatico l’avvio di mongod in fase di accensione
dell’SO:</p>
<p><code>sudo systemctl enable mongod.service</code></p></li>
<li><p>runnare lo script genbank.py, specificando un file .gbk sorgente
e selezionando il tipo di database desiderato (-n per MongoDB, -s per
SQLite3) [DEPRECATED]</p></li>
<li><p>attendere il termine del salvataggio nel database locale</p></li>
<li><p>aprire il server, digitando su console all’interno della cartella
nodeServer: <code>npm install --save express</code> per il primo avvio;
<code>node server.js</code></p></li>
<li><p>dopodiché aprire il browser e cercare il seguente indirizzo:</p>
<p>localhost:5173</p></li>
</ul>
<p>è ora possibile effettuare le query di interesse, ma va inizializzato
il database locale.</p>
<p>Per avere le specie con i dati genomici (versione MongoDB):</p>
<p><code>python3 genbank.py --file "nomefile".gbk -n</code></p>
<p>Per la taxonomy:</p>
<p><code>python3 genbank.py --find -m --email</code> (inserire email per
accesso su ncbi, ammesso che si abbia accesso)</p>
<p>NOTA: a ogni riavvio del sistema, è necessario riattivare mongodb
(<code>systemctl start mongod</code>) e il server
(<code>node server.js</code>)</p>
<h1 id="concetti-utili">Concetti utili:</h1>
<p>MONGODB conserva i dati implicitamente in una memoria virtuale. Per
trasportarli da un sistema a un altro è necessario utilizzare il comando
mongodump per generare una cartella contenente il db di interesse e sudo
mongorestore sul file bson generato dal mongodump una volta importata la
cartella generata.</p>
<h1 id="le-collections-trattate">Le collections trattate:</h1>
<ul>
<li>nucleotide_data contiene tutti i dati delle microalghe che sono
riscontrabili su NCBI alla voce Nucleotide;</li>
<li>nucleotide_basic contiene i dati, allegeriti soltanto a nome e
NCBI_ID, per questioni di performance quando si effettuano query di
conteggio;</li>
<li>taxonomy_data contiene tutti i dati delle microalghe che sono
riscontrabili su NCBI ala voce Taxonomy;</li>
<li>taxonomy_tree contiene i dati parsati per singola specie,
mostrandone di volta in volta i link di Lineage.</li>
</ul>
<h1 id="la-rappresentazione-fasta-e-fastq">La rappresentazione FASTA e
FASTQ</h1>
<p>FASTA conserva soltanto la Sequenza di nucleotidi o amminoacidi,
codificando ogni gene in singole lettere per indice di posizione. Nella
rappresentazione in Genbank, troviamo tale dato nel file JSON che
salviamo in locale, alla voce translation per ogni Coding Sequence sotto
ogni Specie, secondo la seguente gerarchia:</p>
<h1 id="specie">SPECIE</h1>
<h2 id="features">FEATURES</h2>
<h3 id="cds">CDS</h3>
<h4
id="translationlslavgttitlasyhwll">/translation=“LSLAVGTTITLASYHWLL[…]””</h4>
<p>FASTQ è un “quality score” che associa alla sequenza, per ogni indice
di posizione, un valore qualitativo codificato in ASCII. Un esempio a
seguire: <span class="citation" data-cites="SRR64">@SRR64 [...]</span>
Name Sequence CCTCGTCTA[…] DNA Sequence +SRR64[…] Quality address
BBBBBFFFF[…] Quality Score</p>
<h1 id="sequenziamento-genetico">Sequenziamento genetico</h1>
<p>E’ il processo di determinazione dell’ordine dei nucleotidi (Adenina,
Citosina, Guanina e Timina) che costituiscono il frammento di DNA in
analisi. Le tecniche principali di Squenziamento sono Sanger e NGS
(Illumina ne é un esempio).</p>
<h1 id="allineamento-genetico">Allineamento genetico</h1>
<p>E’ il processo di confronto di due o più sequenze di DNA o proteine
per identificare regioni identiche o simili per individuare eventuali
relazioni funzionali, strutturali o filogenetiche. Le tecniche di
allineamento prevedono il confronto globale e locale. Un’applicazione
algoritmica di allineamento locale è data da Smith Waterman.
Un’applicazione algoritmica di allineamento globale é data da
Needleman-Wunsch.</p>
<figure>
<img src="algaeStruct.png" title="Struct of the Project"
alt="Algae Project Struct" />
<figcaption aria-hidden="true">Algae Project Struct</figcaption>
</figure>
<p>Json Files saved:</p>
<p><code>dataStruct.json</code>:</p>
<div class="sourceCode" id="cb1"><pre
class="sourceCode json"><code class="sourceCode json"><span id="cb1-1"><a href="#cb1-1" aria-hidden="true" tabindex="-1"></a>    <span class="fu">{</span></span>
<span id="cb1-2"><a href="#cb1-2" aria-hidden="true" tabindex="-1"></a>        <span class="dt">&quot;struct&quot;</span><span class="fu">:</span><span class="ot">[</span></span>
<span id="cb1-3"><a href="#cb1-3" aria-hidden="true" tabindex="-1"></a>            <span class="st">&quot;Name&quot;</span><span class="er">:</span> <span class="st">&quot;LOCUS&quot;</span><span class="ot">,</span></span>
<span id="cb1-4"><a href="#cb1-4" aria-hidden="true" tabindex="-1"></a>            <span class="st">&quot;Id&quot;</span><span class="er">:</span> <span class="st">&quot;VERSION&quot;</span><span class="ot">,</span></span>
<span id="cb1-5"><a href="#cb1-5" aria-hidden="true" tabindex="-1"></a>            <span class="st">&quot;Seq_Hex&quot;</span><span class="er">:</span> <span class="st">&quot;SHA256HEX&quot;</span><span class="ot">,</span></span>
<span id="cb1-6"><a href="#cb1-6" aria-hidden="true" tabindex="-1"></a>            <span class="st">&quot;Description&quot;</span><span class="er">:</span> <span class="st">&quot;DEFINITION&quot;</span><span class="ot">,</span></span>
<span id="cb1-7"><a href="#cb1-7" aria-hidden="true" tabindex="-1"></a>            <span class="st">&quot;Features&quot;</span><span class="er">:</span> <span class="ot">[</span></span>
<span id="cb1-8"><a href="#cb1-8" aria-hidden="true" tabindex="-1"></a>                <span class="fu">{</span></span>
<span id="cb1-9"><a href="#cb1-9" aria-hidden="true" tabindex="-1"></a>                    <span class="dt">&quot;Type&quot;</span><span class="fu">:</span> <span class="st">&quot;typeFeatures&quot;</span><span class="fu">,</span></span>
<span id="cb1-10"><a href="#cb1-10" aria-hidden="true" tabindex="-1"></a>                    <span class="dt">&quot;Others&quot;</span><span class="fu">:</span> <span class="st">&quot;others info of the feature&quot;</span><span class="fu">,</span></span>
<span id="cb1-11"><a href="#cb1-11" aria-hidden="true" tabindex="-1"></a>                    <span class="dt">&quot;Location&quot;</span><span class="fu">:</span> <span class="ot">[</span></span>
<span id="cb1-12"><a href="#cb1-12" aria-hidden="true" tabindex="-1"></a>                        <span class="dv">0</span><span class="ot">,</span></span>
<span id="cb1-13"><a href="#cb1-13" aria-hidden="true" tabindex="-1"></a>                        <span class="dv">0</span></span>
<span id="cb1-14"><a href="#cb1-14" aria-hidden="true" tabindex="-1"></a>                    <span class="ot">]</span></span>
<span id="cb1-15"><a href="#cb1-15" aria-hidden="true" tabindex="-1"></a>                <span class="fu">}</span></span>
<span id="cb1-16"><a href="#cb1-16" aria-hidden="true" tabindex="-1"></a>            <span class="ot">]</span></span>
<span id="cb1-17"><a href="#cb1-17" aria-hidden="true" tabindex="-1"></a>        <span class="ot">]</span><span class="fu">,</span></span>
<span id="cb1-18"><a href="#cb1-18" aria-hidden="true" tabindex="-1"></a>        <span class="dt">&quot;hex&quot;</span><span class="fu">:</span><span class="ot">[</span></span>
<span id="cb1-19"><a href="#cb1-19" aria-hidden="true" tabindex="-1"></a>            <span class="fu">{</span></span>
<span id="cb1-20"><a href="#cb1-20" aria-hidden="true" tabindex="-1"></a>                <span class="dt">&quot;Seq_Hex&quot;</span><span class="fu">:</span><span class="st">&quot;SHA256HEX&quot;</span><span class="fu">,</span></span>
<span id="cb1-21"><a href="#cb1-21" aria-hidden="true" tabindex="-1"></a>                <span class="dt">&quot;Seq_Raw&quot;</span><span class="fu">:</span><span class="st">&quot;ORIGIN&quot;</span></span>
<span id="cb1-22"><a href="#cb1-22" aria-hidden="true" tabindex="-1"></a>            <span class="fu">}</span></span>
<span id="cb1-23"><a href="#cb1-23" aria-hidden="true" tabindex="-1"></a>        <span class="ot">]</span></span>
<span id="cb1-24"><a href="#cb1-24" aria-hidden="true" tabindex="-1"></a>    <span class="fu">}</span></span>
<span id="cb1-25"><a href="#cb1-25" aria-hidden="true" tabindex="-1"></a>    <span class="er">//example</span></span>
<span id="cb1-26"><a href="#cb1-26" aria-hidden="true" tabindex="-1"></a>    <span class="fu">{</span></span>
<span id="cb1-27"><a href="#cb1-27" aria-hidden="true" tabindex="-1"></a>        <span class="dt">&quot;struct&quot;</span><span class="fu">:</span><span class="ot">[</span></span>
<span id="cb1-28"><a href="#cb1-28" aria-hidden="true" tabindex="-1"></a>            <span class="fu">{</span></span>
<span id="cb1-29"><a href="#cb1-29" aria-hidden="true" tabindex="-1"></a>            <span class="dt">&quot;Name&quot;</span><span class="fu">:</span> <span class="st">&quot;OQ443078&quot;</span><span class="fu">,</span></span>
<span id="cb1-30"><a href="#cb1-30" aria-hidden="true" tabindex="-1"></a>            <span class="dt">&quot;Id&quot;</span><span class="fu">:</span> <span class="st">&quot;OQ443078.1&quot;</span><span class="fu">,</span></span>
<span id="cb1-31"><a href="#cb1-31" aria-hidden="true" tabindex="-1"></a>            <span class="dt">&quot;Seq_Hex&quot;</span><span class="fu">:</span> <span class="st">&quot;de5f363dd78bc0f2a534c55783739a672ce136c67b34767a052d8ecd7e4deb5b&quot;</span><span class="fu">,</span></span>
<span id="cb1-32"><a href="#cb1-32" aria-hidden="true" tabindex="-1"></a>            <span class="dt">&quot;Description&quot;</span><span class="fu">:</span> <span class="st">&quot;Bacillus velezensis strain 3(JS) 16S ribosomal RNA gene, partial sequence&quot;</span><span class="fu">,</span></span>
<span id="cb1-33"><a href="#cb1-33" aria-hidden="true" tabindex="-1"></a>            <span class="dt">&quot;Features&quot;</span><span class="fu">:</span> <span class="ot">[</span></span>
<span id="cb1-34"><a href="#cb1-34" aria-hidden="true" tabindex="-1"></a>                <span class="fu">{</span></span>
<span id="cb1-35"><a href="#cb1-35" aria-hidden="true" tabindex="-1"></a>                    <span class="dt">&quot;Type&quot;</span><span class="fu">:</span> <span class="st">&quot;source&quot;</span><span class="fu">,</span></span>
<span id="cb1-36"><a href="#cb1-36" aria-hidden="true" tabindex="-1"></a>                    <span class="dt">&quot;organism&quot;</span><span class="fu">:</span> <span class="st">&quot;Bacillus velezensis&quot;</span><span class="fu">,</span></span>
<span id="cb1-37"><a href="#cb1-37" aria-hidden="true" tabindex="-1"></a>                    <span class="dt">&quot;mol_type&quot;</span><span class="fu">:</span> <span class="st">&quot;genomic DNA&quot;</span><span class="fu">,</span></span>
<span id="cb1-38"><a href="#cb1-38" aria-hidden="true" tabindex="-1"></a>                    <span class="dt">&quot;strain&quot;</span><span class="fu">:</span> <span class="st">&quot;3(JS)&quot;</span><span class="fu">,</span></span>
<span id="cb1-39"><a href="#cb1-39" aria-hidden="true" tabindex="-1"></a>                    <span class="dt">&quot;db_xref&quot;</span><span class="fu">:</span> <span class="st">&quot;taxon:492670&quot;</span><span class="fu">,</span></span>
<span id="cb1-40"><a href="#cb1-40" aria-hidden="true" tabindex="-1"></a>                    <span class="dt">&quot;country&quot;</span><span class="fu">:</span> <span class="st">&quot;India: Vikramgad&quot;</span><span class="fu">,</span></span>
<span id="cb1-41"><a href="#cb1-41" aria-hidden="true" tabindex="-1"></a>                    <span class="dt">&quot;lat_lon&quot;</span><span class="fu">:</span> <span class="st">&quot;19.80 N 73.09 E&quot;</span><span class="fu">,</span></span>
<span id="cb1-42"><a href="#cb1-42" aria-hidden="true" tabindex="-1"></a>                    <span class="dt">&quot;collected_by&quot;</span><span class="fu">:</span> <span class="st">&quot;Vir Acharya&quot;</span><span class="fu">,</span></span>
<span id="cb1-43"><a href="#cb1-43" aria-hidden="true" tabindex="-1"></a>                    <span class="dt">&quot;identified_by&quot;</span><span class="fu">:</span> <span class="st">&quot;Raunak Giri&quot;</span><span class="fu">,</span></span>
<span id="cb1-44"><a href="#cb1-44" aria-hidden="true" tabindex="-1"></a>                    <span class="dt">&quot;Location&quot;</span><span class="fu">:</span> <span class="ot">[</span></span>
<span id="cb1-45"><a href="#cb1-45" aria-hidden="true" tabindex="-1"></a>                        <span class="dv">0</span><span class="ot">,</span></span>
<span id="cb1-46"><a href="#cb1-46" aria-hidden="true" tabindex="-1"></a>                        <span class="dv">936</span></span>
<span id="cb1-47"><a href="#cb1-47" aria-hidden="true" tabindex="-1"></a>                    <span class="ot">]</span></span>
<span id="cb1-48"><a href="#cb1-48" aria-hidden="true" tabindex="-1"></a>                <span class="fu">}</span><span class="ot">,</span></span>
<span id="cb1-49"><a href="#cb1-49" aria-hidden="true" tabindex="-1"></a>                <span class="fu">{</span></span>
<span id="cb1-50"><a href="#cb1-50" aria-hidden="true" tabindex="-1"></a>                    <span class="dt">&quot;Type&quot;</span><span class="fu">:</span> <span class="st">&quot;rRNA&quot;</span><span class="fu">,</span></span>
<span id="cb1-51"><a href="#cb1-51" aria-hidden="true" tabindex="-1"></a>                    <span class="dt">&quot;product&quot;</span><span class="fu">:</span> <span class="st">&quot;16S ribosomal RNA&quot;</span><span class="fu">,</span></span>
<span id="cb1-52"><a href="#cb1-52" aria-hidden="true" tabindex="-1"></a>                    <span class="dt">&quot;Location&quot;</span><span class="fu">:</span> <span class="ot">[</span></span>
<span id="cb1-53"><a href="#cb1-53" aria-hidden="true" tabindex="-1"></a>                        <span class="dv">0</span><span class="ot">,</span></span>
<span id="cb1-54"><a href="#cb1-54" aria-hidden="true" tabindex="-1"></a>                        <span class="dv">936</span></span>
<span id="cb1-55"><a href="#cb1-55" aria-hidden="true" tabindex="-1"></a>                    <span class="ot">]</span></span>
<span id="cb1-56"><a href="#cb1-56" aria-hidden="true" tabindex="-1"></a>                <span class="fu">}</span></span>
<span id="cb1-57"><a href="#cb1-57" aria-hidden="true" tabindex="-1"></a>            <span class="ot">]</span></span>
<span id="cb1-58"><a href="#cb1-58" aria-hidden="true" tabindex="-1"></a>        <span class="fu">}</span></span>
<span id="cb1-59"><a href="#cb1-59" aria-hidden="true" tabindex="-1"></a>        <span class="ot">]</span><span class="fu">,</span></span>
<span id="cb1-60"><a href="#cb1-60" aria-hidden="true" tabindex="-1"></a>        <span class="dt">&quot;hex&quot;</span><span class="fu">:</span><span class="ot">[</span></span>
<span id="cb1-61"><a href="#cb1-61" aria-hidden="true" tabindex="-1"></a>            <span class="fu">{</span></span>
<span id="cb1-62"><a href="#cb1-62" aria-hidden="true" tabindex="-1"></a>                <span class="dt">&quot;Seq_Hex&quot;</span><span class="fu">:</span><span class="st">&quot;de5f363dd78bc0f2a534c55783739a672ce136c67b34767a052d8ecd7e4deb5b&quot;</span><span class="fu">,</span></span>
<span id="cb1-63"><a href="#cb1-63" aria-hidden="true" tabindex="-1"></a>                <span class="dt">&quot;Seq_Raw&quot;</span><span class="fu">:</span><span class="st">&quot;ACCTGCCTGTAAGACTGGGATAACTCCGGGAAACCGGGGCTAATACCGGATGGTTGTCTGAACCGCATGGTTCAGACATAAAAGGTGGCTTCGGCTACCACTTACAGATGGACCCGCGGCGCATTAGCTAGTTGGTGAGGTAACGGCTCACCAAGGCGACGATGCGTAGCCGACCTGAGAGGGTGATCGGCCACACTGGGACTGAGACACGGCCCAGACTCCTACGGGAGGCAGCAGTAGGGAATCTTCCGCAATGGACGAAAGTCTGACGGAGCAACGCCGCGTGAGTGATGAAGGTTTTCGGATCGTAAAGCTCTGTTGTTAGGGAAGAACAAGTGCCGTTCAAATAGGGCGGCACCTTGACGGTACCTAACCAGAAAGCCACGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGTGGCAAGCGTTGTCCGGAATTATTGGGCGTAAAGGGCTCGCAGGCGGTTTCTTAAGTCTGATGTGAAAGCCCCCGGCTCAACCGGGGAGGGTCATTGGAAACTGGGGAACTTGAGTGCAGAAGAGGAGAGTGGAATTCCACGTGTAGCGGTGAAATGCGTAGAGATGTGGAGGAACACCAGTGGCGAAGGCGACTCTCTGGTCTGTAACTGACGCTGAGGAGCGAAAGCGTGGGGAGCGAACAGGATTAGATACCCTGGTAGTCCACGCCGTAAACGATGAGTGCTAAGTGTTAGGGGGTTTCCGCCCCTTAGTGCTGCAGCTAACGCATTAAGCACTCCGCCTGGGGAGTACGGTCGCAAGACTGAAACTCAAAGGAATTGACGGGGGCCCGCACAAGCGGTGGAGCATGTGGTTTAATTCGAAGCAACGCGAAGAACCTTACCAGGTCTTGACATCCTCTGACAATCCTAGAGATAGGACGTCCCCTTCGGGGGCAGAGTGACAGGTGGTGC&quot;</span></span>
<span id="cb1-64"><a href="#cb1-64" aria-hidden="true" tabindex="-1"></a>            <span class="fu">}</span></span>
<span id="cb1-65"><a href="#cb1-65" aria-hidden="true" tabindex="-1"></a>        <span class="ot">]</span></span>
<span id="cb1-66"><a href="#cb1-66" aria-hidden="true" tabindex="-1"></a>    <span class="fu">}</span></span></code></pre></div>
<p><code>dataSource.json</code></p>
<div class="sourceCode" id="cb2"><pre
class="sourceCode json"><code class="sourceCode json"><span id="cb2-1"><a href="#cb2-1" aria-hidden="true" tabindex="-1"></a>    <span class="fu">{</span></span>
<span id="cb2-2"><a href="#cb2-2" aria-hidden="true" tabindex="-1"></a>        <span class="dt">&quot;OrganismName&quot;</span><span class="fu">:</span><span class="ot">[</span></span>
<span id="cb2-3"><a href="#cb2-3" aria-hidden="true" tabindex="-1"></a>            <span class="st">&quot;VERSION&quot;</span></span>
<span id="cb2-4"><a href="#cb2-4" aria-hidden="true" tabindex="-1"></a>        <span class="ot">]</span></span>
<span id="cb2-5"><a href="#cb2-5" aria-hidden="true" tabindex="-1"></a>    <span class="fu">}</span></span>
<span id="cb2-6"><a href="#cb2-6" aria-hidden="true" tabindex="-1"></a>    <span class="er">//example</span></span>
<span id="cb2-7"><a href="#cb2-7" aria-hidden="true" tabindex="-1"></a>    <span class="fu">{</span></span>
<span id="cb2-8"><a href="#cb2-8" aria-hidden="true" tabindex="-1"></a>        <span class="dt">&quot;Danio rerio&quot;</span><span class="fu">:</span> <span class="ot">[</span></span>
<span id="cb2-9"><a href="#cb2-9" aria-hidden="true" tabindex="-1"></a>            <span class="st">&quot;NM_001327985.1&quot;</span><span class="ot">,</span></span>
<span id="cb2-10"><a href="#cb2-10" aria-hidden="true" tabindex="-1"></a>            <span class="st">&quot;NM_212741.1&quot;</span><span class="ot">,</span></span>
<span id="cb2-11"><a href="#cb2-11" aria-hidden="true" tabindex="-1"></a>            <span class="st">&quot;NM_001128675.2&quot;</span><span class="ot">,</span></span>
<span id="cb2-12"><a href="#cb2-12" aria-hidden="true" tabindex="-1"></a>            <span class="st">&quot;NM_001081554.1&quot;</span><span class="ot">,</span></span>
<span id="cb2-13"><a href="#cb2-13" aria-hidden="true" tabindex="-1"></a>            <span class="st">&quot;NM_199618.1&quot;</span><span class="ot">,</span></span>
<span id="cb2-14"><a href="#cb2-14" aria-hidden="true" tabindex="-1"></a>            <span class="st">&quot;NM_001303262.1&quot;</span><span class="ot">,</span></span>
<span id="cb2-15"><a href="#cb2-15" aria-hidden="true" tabindex="-1"></a>            <span class="st">&quot;NM_001137555.3&quot;</span><span class="ot">,</span></span>
<span id="cb2-16"><a href="#cb2-16" aria-hidden="true" tabindex="-1"></a>            <span class="st">&quot;NM_001020798.1&quot;</span></span>
<span id="cb2-17"><a href="#cb2-17" aria-hidden="true" tabindex="-1"></a>        <span class="ot">]</span></span>
<span id="cb2-18"><a href="#cb2-18" aria-hidden="true" tabindex="-1"></a>    <span class="fu">}</span></span></code></pre></div>
<p><code>algaeSource.json</code></p>
<div class="sourceCode" id="cb3"><pre
class="sourceCode json"><code class="sourceCode json"><span id="cb3-1"><a href="#cb3-1" aria-hidden="true" tabindex="-1"></a>    <span class="fu">{</span></span>
<span id="cb3-2"><a href="#cb3-2" aria-hidden="true" tabindex="-1"></a>        <span class="dt">&quot;OrganismName&quot;</span><span class="fu">:</span><span class="ot">[</span></span>
<span id="cb3-3"><a href="#cb3-3" aria-hidden="true" tabindex="-1"></a>            <span class="st">&quot;VERSION&quot;</span></span>
<span id="cb3-4"><a href="#cb3-4" aria-hidden="true" tabindex="-1"></a>        <span class="ot">]</span></span>
<span id="cb3-5"><a href="#cb3-5" aria-hidden="true" tabindex="-1"></a>    <span class="fu">}</span></span>
<span id="cb3-6"><a href="#cb3-6" aria-hidden="true" tabindex="-1"></a>    <span class="er">//example</span></span>
<span id="cb3-7"><a href="#cb3-7" aria-hidden="true" tabindex="-1"></a>    <span class="fu">{</span></span>
<span id="cb3-8"><a href="#cb3-8" aria-hidden="true" tabindex="-1"></a>        <span class="dt">&quot;Bracteacoccus aggregatus&quot;</span><span class="fu">:</span> <span class="ot">[</span></span>
<span id="cb3-9"><a href="#cb3-9" aria-hidden="true" tabindex="-1"></a>            <span class="st">&quot;MZ090013.1&quot;</span><span class="ot">,</span></span>
<span id="cb3-10"><a href="#cb3-10" aria-hidden="true" tabindex="-1"></a>            <span class="st">&quot;MZ067570.1&quot;</span><span class="ot">,</span></span>
<span id="cb3-11"><a href="#cb3-11" aria-hidden="true" tabindex="-1"></a>            <span class="st">&quot;MH205944.1&quot;</span><span class="ot">,</span></span>
<span id="cb3-12"><a href="#cb3-12" aria-hidden="true" tabindex="-1"></a>            <span class="st">&quot;MH703758.1&quot;</span><span class="ot">,</span></span>
<span id="cb3-13"><a href="#cb3-13" aria-hidden="true" tabindex="-1"></a>            <span class="st">&quot;MH703740.1&quot;</span></span>
<span id="cb3-14"><a href="#cb3-14" aria-hidden="true" tabindex="-1"></a>        <span class="ot">]</span></span>
<span id="cb3-15"><a href="#cb3-15" aria-hidden="true" tabindex="-1"></a>    <span class="fu">}</span></span></code></pre></div>
<p><code>microAlgaeSource.json</code></p>
<div class="sourceCode" id="cb4"><pre
class="sourceCode json"><code class="sourceCode json"><span id="cb4-1"><a href="#cb4-1" aria-hidden="true" tabindex="-1"></a>    <span class="fu">{</span></span>
<span id="cb4-2"><a href="#cb4-2" aria-hidden="true" tabindex="-1"></a>        <span class="dt">&quot;OrganismName&quot;</span><span class="fu">:</span><span class="ot">[</span></span>
<span id="cb4-3"><a href="#cb4-3" aria-hidden="true" tabindex="-1"></a>            <span class="st">&quot;VERSION&quot;</span></span>
<span id="cb4-4"><a href="#cb4-4" aria-hidden="true" tabindex="-1"></a>        <span class="ot">]</span></span>
<span id="cb4-5"><a href="#cb4-5" aria-hidden="true" tabindex="-1"></a>    <span class="fu">}</span></span>
<span id="cb4-6"><a href="#cb4-6" aria-hidden="true" tabindex="-1"></a>    <span class="er">//example</span></span>
<span id="cb4-7"><a href="#cb4-7" aria-hidden="true" tabindex="-1"></a>    <span class="fu">{</span></span>
<span id="cb4-8"><a href="#cb4-8" aria-hidden="true" tabindex="-1"></a>        <span class="dt">&quot;Ankistrodesmus fusiformis&quot;</span><span class="fu">:</span> <span class="ot">[</span></span>
<span id="cb4-9"><a href="#cb4-9" aria-hidden="true" tabindex="-1"></a>            <span class="st">&quot;OM683277.1&quot;</span><span class="ot">,</span></span>
<span id="cb4-10"><a href="#cb4-10" aria-hidden="true" tabindex="-1"></a>            <span class="st">&quot;OM683275.1&quot;</span></span>
<span id="cb4-11"><a href="#cb4-11" aria-hidden="true" tabindex="-1"></a>        <span class="ot">]</span></span>
<span id="cb4-12"><a href="#cb4-12" aria-hidden="true" tabindex="-1"></a>    <span class="fu">}</span></span></code></pre></div>

		<ErrorMsg v-if="errormsg" :msg="errormsg"></ErrorMsg>
	</div>
</template>
