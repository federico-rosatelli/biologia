const express = require('express')
const app = express()
const path = require('path');

const methodOverride = require('method-override')
const bodyParser = require('body-parser')
const http = require('http');

const server = http.createServer(app);

app.use(express.static('public'))
app.use(express.static(path.join(__dirname, 'views/')));

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




async function finder(query) {
  const client = await MongoClient.connect(uri, { useNewUrlParser: true ,useUnifiedTopology: true })
  .catch(err => { console.log(err); });
  try {
    const database = client.db('Biologia');
    const user = database.collection('taxonomy_data');
    const result = await user.find(query).toArray();
    await client.close();
    return result;
  }
  catch{
    console.log("error qui");
  }
}

app.get('/',async(req,res)=>{
    let query = ""
    let rank = ""
    if (req.query.search == ""){
      let info = {Lineage:{$regex:`${query}`}};
      var find = await finder(info)
        let count = {
          superkingdom:{},
          kingdom:{},
          phylum:{},
          clade:{},
          class:{},
          clade:{},
          order:{},
          family:{},
          genus:{},
          "no rank":{},
          subfamily:{},
          subclass:{},
          species:{}
        }
        for (let i = 0; i < find.length; i++) {
          let data = find[i]
          for (let j = 0; j < data.LineageEx.length; j++) {
            if(data.LineageEx[j].ScientificName in count[data.LineageEx[j].Rank]){
              count[data.LineageEx[j].Rank][data.LineageEx[j].ScientificName] += 1;
            }
            else{
              count[data.LineageEx[j].Rank][data.LineageEx[j].ScientificName] = 1;
            }
          }
        }
        for (let i = 0; i < find.length; i++) {
          let data = find[i]
          for (let j = 0; j < data.LineageEx.length; j++) {

              find[i].LineageEx[j]["tot"] = count[data.LineageEx[j].Rank][data.LineageEx[j].ScientificName]
          }
        }
        res.render('home',{find})
    }
    else if (req.query.search) {
        
        query = req.query.search;
        let info = {}
        if (req.query.rank){
          rank = req.query.rank;

          info = {Lineage:{$regex:`${query}`},LineageEx:{$elemMatch:{Rank:rank,ScientificName:query}}};
        }
        else{
          info = {Lineage:{$regex:`${query}`}}
        }
        var find = await finder(info)
        let count = {
          superkingdom:{},
          kingdom:{},
          phylum:{},
          clade:{},
          class:{},
          clade:{},
          order:{},
          family:{},
          genus:{},
          "no rank":{},
          subfamily:{},
          subclass:{}
        }
        for (let i = 0; i < find.length; i++) {
          let data = find[i]
          for (let j = 0; j < data.LineageEx.length; j++) {
            if(data.LineageEx[j].ScientificName in count[data.LineageEx[j].Rank]){
              count[data.LineageEx[j].Rank][data.LineageEx[j].ScientificName] += 1;
            }
            else{
              count[data.LineageEx[j].Rank][data.LineageEx[j].ScientificName] = 1;
            }
          }
        }
        for (let i = 0; i < find.length; i++) {
          let data = find[i]
          for (let j = 0; j < data.LineageEx.length; j++) {

              find[i].LineageEx[j]["tot"] = count[data.LineageEx[j].Rank][data.LineageEx[j].ScientificName]
          }
        }
        res.render('home',{find})
    }
    else{
        find = null
        res.render('home',{find:null}) 
    }


})


server.listen(process.env.PORT || 3000)

