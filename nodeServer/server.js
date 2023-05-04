const express = require('express')
const app = express()
const path = require('path');

const methodOverride = require('method-override')
const bodyParser = require('body-parser')
const http = require('http');

const server = http.createServer(app);

app.use(express.static('public'))
app.use(express.static(path.join(__dirname, 'views/')));  // Permette al codice CSS di essere applicato alla pagina web.

app.use(bodyParser.urlencoded({
  extended: true
}));
app.use(bodyParser.json());

app.set('view engine','ejs')
app.use(express.urlencoded({
  extended: false
}))

app.use(methodOverride('_method'))


const { MongoClient } = require('mongodb');

const uri = "mongodb://localhost:27017";


// Finder methods

async function finder(query, collection) {
  console.log('Begin of search in Biologia DB, in collection '+collection)
  const client = await MongoClient.connect(uri, { useNewUrlParser: true ,useUnifiedTopology: true })
  .catch(err => { console.log(err); });
  try {
    const database = client.db('Biologia');
    const user = database.collection(collection);
    const result = await user.find(query).toArray();
    await client.close();
    //console.log('Return result '+ result + "\n")
    return result;
  }
  catch(e){
    console.log(`exception in finder(query): Exception ${e}`);
    return [];
  }
}

async function finderOne(query, collection) {
  const client = await MongoClient.connect(uri, { useNewUrlParser: true ,useUnifiedTopology: true })
  .catch(err => { console.log(err); });
  try {
    const database = client.db('Biologia');
    const user = database.collection(collection);
    const result = await user.findOne(query);
    await client.close();
    return result;
  }
  catch{
    console.log("exception in finder(query)");
  }
}

app.get('/',async(req, res) => {
    let query = ""
    //let taxonomy = ""
    let data = ""
    let nucleotide = null
    let taxonomy = null
    // let rank = ""
    if (req.query.data) {
      console.log("Data in input. Begin process...")
      data = req.query.data
      let info = {Lineage:{$regex:`^${data}`,$options : 'i'}}; // potrebbe essere cosi? boh proviamo
      var find1 = await finder(info, 'taxonomy_data')
      var array = []
      for (let i = 0; i < find1.length; i++) {
        let info2 = {Features:{$elemMatch:{Type:'source', db_xref:"taxon:"+find1[i].TaxId}}}
        // qua si cerca un nucleotide in particolare giusto?
        //io stavo pensando ad una cosa del tipo
        // Chlorella Vulgaris: [
        //    dato di nucleotide 1: []
        //    dato di nucleotide 2: []
        // ]
        // è da fare ma per ora stavamo facendo un array e basta perché stiamo morendo e volevamo fare in fretta
        // finder trova una lista di elementi mentre se ce ne serve solo uno dobbiamo utilizzare un metodo diverso
        // ho creato la funzione finderOne che trova una solo occorrenza invece di many
        // ahahha ok


        var find2 = await finder(info2, 'nucleotide_data')
        if (find2[i] != undefined) {
          for (let o = 0; o < find2.length; o++) {
          // FOR DEBUGGING PURPOSE
          // if (find2[i] != undefined) {
          //   console.log(find2[i].Id)
          // }
          // else {
          //   console.log("Found undefined! "+o)
          // }
            array.push(find2[i])
          }
        }        
      }
      console.log("Query ended. Result displayed at localhost:3000")
      res.render('home', {find:taxonomy,nucleotide:array})    // {find:array} prende array e gli dice che si chiama find  
    }
    if (req.query.taxonomy == ''){
      console.log("No values entered. Begin process...")
      let info = {Lineage:{$regex:`${query}`}};
      var find = await finder(info, 'taxonomy_data')
        // let count = {
        //   superkingdom:{},
        //   kingdom:{},
        //   phylum:{},
        //   clade:{},
        //   class:{},
        //   clade:{},
        //   order:{},
        //   family:{},
        //   genus:{},
        //   "no rank":{},
        //   subfamily:{},
        //   subclass:{},
        //   species:{}
        // }
        // for (let i = 0; i < find.length; i++) {
        //   let data = find[i]
        //   for (let j = 0; j < data.LineageEx.length; j++) {
        //     if(data.LineageEx[j].ScientificName in count[data.LineageEx[j].Rank]){
        //       count[data.LineageEx[j].Rank][data.LineageEx[j].ScientificName] += 1;
        //     }
        //     else{
        //       count[data.LineageEx[j].Rank][data.LineageEx[j].ScientificName] = 1;
        //     }
        //   }
        // }
        // for (let i = 0; i < find.length; i++) {
        //   let data = find[i]
        //   for (let j = 0; j < data.LineageEx.length; j++) {

        //       find[i].LineageEx[j]["tot"] = count[data.LineageEx[j].Rank][data.LineageEx[j].ScientificName]
        //   }
        // }
        console.log("Query ended. Result displayed at localhost:3000")
        res.render('home',{find,nucleotide})
    }
    else if (req.query.taxonomy) {
        console.log("Detected some key values in Search Field Taxonomy. Begin process...")
        query = req.query.taxonomy;
        let info = {}

        //
        //  TO DO: dare in input al DB un testo e tradurlo in lowercase
        //  Da aggiungere la ricerca per specie
        //
        if (req.query.rank){
          console.log("Rank in input.")
          rank = req.query.rank;
          info = {LineageEx:{$elemMatch:{Rank:rank,ScientificName:query}}};
        }
        else{
          console.log("No Rank in input.")
          info = {Lineage:{$regex:`${query}`}}
        }
        var find = await finder(info, 'taxonomy_data')
        // let count = {
        //   superkingdom:{},
        //   kingdom:{},
        //   phylum:{},
        //   clade:{},
        //   class:{},
        //   clade:{},
        //   order:{},
        //   family:{},
        //   genus:{},
        //   "no rank":{},
        //   subfamily:{},
        //   subclass:{}
        // }
        // for (let i = 0; i < find.length; i++) {
        //   let data = find[i]
        //   for (let j = 0; j < data.LineageEx.length; j++) {
        //     if(data.LineageEx[j].ScientificName in count[data.LineageEx[j].Rank]){
        //       count[data.LineageEx[j].Rank][data.LineageEx[j].ScientificName] += 1;
        //     }
        //     else{
        //       count[data.LineageEx[j].Rank][data.LineageEx[j].ScientificName] = 1;
        //     }
        //   }
        // }
        // for (let i = 0; i < find.length; i++) {
        //   let data = find[i]
        //   for (let j = 0; j < data.LineageEx.length; j++) {

        //       find[i].LineageEx[j]["tot"] = count[data.LineageEx[j].Rank][data.LineageEx[j].ScientificName]
        //   }
        // }
        console.log("Query ended. Result displayed at localhost:3000")
        res.render('home',{find,nucleotide})
    }
    else if (req.query.sequenceID) {
      console.log("Detected some key values in Search Field Product. Begin process...")
      query = req.query.sequenceID;
      // if (query == ""){
      //   query = "blabla";
      // }
      let type = req.query.id;
      console.log(type);
      let id = null
      if (id == "genome"){
        id = "locus_tag";
      }
      //let type = req.query.id
      //let info = {Features:{$elemMatch:{Type:'source', db_xref:"taxon:"+find1[i].TaxId}}}
      // query = "small subunit ribosomal RNA"
      let dataType = {
      }
      dataType[type] = {$regex:`${query}`};

      if (type == "genome"){
        dataType["Type"] = "CDS";
        dataType["$exists"] = {"locus_tag":1}
      }

      let info = {Features:{$elemMatch:dataType}};
      // let info = {"Features.$.product":`${query}`}
      nucleotide = await finder(info, 'nucleotide_data');
      
      res.render('home', {find,nucleotide});
    }
    else{
        find = null
        console.log("Result empty. Displayed at localhost:3000")
        res.render('home',{find:null,nucleotide})
    }
    console.log("Return at Home\n")
    

})


server.listen(process.env.PORT || 3000)

