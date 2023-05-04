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
    if (result == null) {
      console.log('Value not found\n');
    }
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
      let info = {ScientificName:{$regex:`${data}`}}; // potrebbe essere cosi? boh proviamo ,$options : 'i'
      var find1 = await finder(info, 'taxonomy_data')
      var array = []
      for (let i = 0; i < find1.length; i++) {
        let info2 = {Features:{$elemMatch:{Type:'source', db_xref:"taxon:"+find1[i].TaxId}}}
        var find2 = await finder(info2, 'nucleotide_data')
        let ids = []
        let array = []
        if (!info2){
          res.render('home', {find:taxonomy,nucleotide,error:true})
          return
        }
        for (let k = 0; k < find2.length; k++) {
          if (!ids.includes(find2[k].Id)){
            try {
              find2[i]["Lineage"] = find1[i].Lineage
              ids.push(find2[k].Id)
              array.push(find2[k])
            }
            catch {
              console.log("Critical Error in search.")
            }
            
          }
        }
        // if (find2[i] != undefined) {
        //   for (let o = 0; o < find2.length; o++) {
        //   // FOR DEBUGGING PURPOSE
        //   // if (find2[i] != undefined) {
        //   //   console.log(find2[i].Id)
        //   // }
        //   // else {
        //   //   console.log("Found undefined! "+o)
        //   // }
        //     array.push(find2[i])
        //   }
        // }        
      }
      //console.log(array);
      console.log("Query ended. Result displayed at localhost:3000")
      res.render('home', {find:taxonomy,nucleotide:find2,error:false})    // {find:array} prende array e gli dice che si chiama find  
      return
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
        res.render('home',{find,nucleotide,error:false})
        return
    }
    //  TO DO: Parsing di Taxonomy incompleto, fixare genbank.py per eliminare errore HTTP ERROR 400 quando trova undefined
    //  Da aggiungere la ricerca per specie
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
        res.render('home',{find,nucleotide,error:false})
        return
    }
    else if (req.query.sequenceID) {
      console.log("Detected some key values in Search Field Product. Begin process...")
      query = req.query.sequenceID;
      if (query == ""){
        query = "blabla";
      }
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
      
      if (type == "any"){
        dataType["$or"] = []
        info1 = {}
        info1["Type"] = "CDS";
        info1["locus_tag"] = {$regex:`${query}`}
        info2 = {}
        info2["protein_id"] = {$regex:`${query}`};
        info3 = {}
        info3["product"] = {$regex:`${query}`};
        dataType["$or"].push(info1,info2,info3)
      }

      else if (type == "genome"){
        dataType["Type"] = "CDS";
        dataType["locus_tag"] = {$regex:`${query}`}
      }
      else{
        dataType[type] = {$regex:`${query}`};
      }

      console.log(dataType);

      let info = {Features:{$elemMatch:dataType}};
      // let info = {"Features.$.product":`${query}`}
      nucleotide = await finder(info, 'nucleotide_data');
      results = []
      let organismCount = {}
      let proteinCount = {}
      for (let i = 0; i < nucleotide.length; i++) {
        let res = nucleotide[i]
        for (let k = 0; k < res.Features.length; k++) {
          let feature = res.Features[k]
          let dataPush = {
            
          }
          if (feature.Type == "source"){
            dataPush.organism = feature.organism
            if (!organismCount[feature.organism]){
              organismCount[feature.organism] = 1
            }
          }
          else if (feature.Type == "CDS"){
            dataPush.protein = feature.protein_id
          }
          dataPush.organismCount = 0
          dataPush.proteinCount = 0
          results.push(dataPush)
        }
        
      }
      res.render('home', {find,nucleotide:results,error:false});
      return
    }
    else if(req.query.sequenceID == "" || req.query.data == ""){
        find = null
        console.log("Result empty. Displayed at localhost:3000")
        res.render('home',{find:null,nucleotide,error:true})
        return
    }
    console.log("Return at Home\n")
    res.render('home',{find:null,nucleotide:null,error:false})
    return

})


server.listen(process.env.PORT || 3000)

