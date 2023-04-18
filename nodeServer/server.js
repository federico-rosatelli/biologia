const express = require('express')
const app = express()
const path = require('path');

const methodOverride = require('method-override')
const bodyParser = require('body-parser')
const http = require('http');

const server = http.createServer(app);

app.use(express.static('public'))

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
    if (req.query.search) {
        query = req.query.search;
        var find = await finder({Lineage:{$regex:`${query}`}})
        
        res.render('home',{find})
    }
    else{
        find = null
        res.render('home',{find:null}) 
    }


})


server.listen(process.env.PORT || 3000)

