const express = require('express')
const app = express()
const path = require('path');
const fs = require('fs')

const methodOverride = require('method-override')
const bodyParser = require('body-parser')
const http = require('http');

const server = http.createServer(app);

app.use(express.static('public'))
app.use(express.static(path.join(__dirname, 'views/')));  // Permette al codice CSS di essere applicato alla pagina web.

app.use(bodyParser.urlencoded({
  extended: true
}));
app.use(express.json());
app.use(express.urlencoded({ extended: true }));
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

app.post('/save',async(req,res) =>{
    let datas = req.body;
    let data = JSON.stringify(datas);
    fs.writeFileSync('public/data.json', data);
    res.setHeader('Content-Length', 0);
    res.download(path.join(__dirname,'public/data.json'));
})

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
              find2[k]["Lineage"] = find1[i].Lineage
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
        let taxon = req.query.taxonomy;

        //
        //  TO DO: dare in input al DB un testo e tradurlo in lowercase
        //  Da aggiungere la ricerca per specie
        //
        // if (req.query.rank){
        //   console.log("Rank in input.")
        //   rank = req.query.rank;
        //   info = {LineageEx:{$elemMatch:{Rank:rank,ScientificName:query}}};
        // }
        // else{
        //   console.log("No Rank in input.")
        //   info = {Lineage:{$regex:`${query}`}}
        // }
        // var find = await finder(info, 'taxonomy_data')
        let filter = req.query.filter
        let search = {}
        if (filter == "scientific_name"){
          search["ScientificName"] = `${taxon}`
        }
        else if (filter == "taxid"){
          search["TaxId"] = `${taxon}`
        }
        else{
          res.render('home',{find:null,nucleotide,error:false})
          return
        }
        let find = await finderOne(search, 'taxonomy_tree')
        if (!find){
          res.render('home',{find:null,nucleotide,error:false})
          return
        }
        let dataReturn = {}
        dataReturn = {
          'TaxId':find.TaxId,
          'ScientificName': find.ScientificName,
          'Rank':find.Rank,
          'SubClasses':[]
        }
        let i1 = 0;
        let i2 = 0;
        let i3 = 0;
        while(i1 <= find.SubClasses.length){
          let dataClass1 = {}
          let tx1 = find.SubClasses[i1]
          if (!tx1){
            i1++;
            continue
          }
          dataClass1 = {
            'TaxId':tx1.TaxId,
            'ScientificName': tx1.ScientificName,
            'Rank':tx1.Rank,
            'SubClasses':[]
          }
          dataReturn.SubClasses.push(dataClass1)
          let filter1 = {TaxId:`${tx1.TaxId}`}
          let find1 = await finderOne(filter1, 'taxonomy_tree')
          if (!find1){
            i1++;
            continue
          }
          while(i2 <= find1.SubClasses.length){
            let dataClass2 = {}
            let tx2 = find1.SubClasses[i2]
            if (!tx2){
              i2++;
              continue
            }
            dataClass2 = {
              'TaxId':tx2.TaxId,
              'ScientificName': tx2.ScientificName,
              'Rank':tx2.Rank,
              'SubClasses':[]
            }
            dataClass1.SubClasses.push(dataClass2)
            let filter2 = {TaxId:`${tx2.TaxId}`}
            let find2 = await finderOne(filter2, 'taxonomy_tree')
            if (!find2){
              i2++;
              continue
            }
            while(i3 <= find2.SubClasses.length){
              let dataClass3 = {}
              let tx3 = find2.SubClasses[i3]
              if (!tx3){
                i3++;
                continue
              }
              dataClass3 = {
                'TaxId':tx3.TaxId,
                'ScientificName': tx3.ScientificName,
                'Rank':tx3.Rank,
              }
              dataClass2.SubClasses.push(dataClass3)
              let filter3 = {TaxId:`${tx3.TaxId}`}
              let find3 = await finderOne(filter3, 'taxonomy_tree')
              if (!find3){
                i3++;
                continue
              }
              
              i3++;
            }
            i3 = 0;
            
            i2++;
          }
          i2 = 0;
          
          i1++;
        }

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
        res.render('home',{find:dataReturn,nucleotide,error:false})
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
      for (let i = 0; i < nucleotide.length; i++) {
        let res = nucleotide[i]
        let dataPush = {}
        dataPush.proteins = []
        for (let k = 0; k < res.Features.length; k++) {
          let feature = res.Features[k]
          
          if (feature.Type == "source"){
            dataPush.organism = feature.organism
            if (!organismCount[feature.organism]){
              organismCount[feature.organism] = 1
            }
            else{
              organismCount[feature.organism] += 1
            }
          }
          else if (feature.Type == "CDS"){
            dataPush.proteins.push(feature.protein_id)
            // if (!proteinCount[feature.protein_id]){
            //   console.log("ENTRA IN CDS");
            //   proteinCount[feature.protein_id] = 1
            // }
            // else{
            //   proteinCount[feature.protein_id] += 1
            // }
          }
          dataPush.organismCount = 0

          results.push(dataPush)
        }
        for (let i = 0; i < results.length; i++) {
          results[i].organismCount = organismCount[results[i].organism]
          results[i].proteinCount = results[i].proteins.length
          
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

