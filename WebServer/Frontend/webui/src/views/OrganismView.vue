<script>
import $ from 'jquery'
export default {
	data: function() {
		return {
            taxId: null,
            typeof: null,
			errormsg: null,
			loading: false,
            response: [],
            nucleo:null,
            dataArray: [],
            otherRouterLink:null,
            taxonomy:null

		}
	},
    
	methods: {
        async refresh() {
            try{
                let response = await this.$axios.get(`/organism/${this.taxId}/${this.typeof}`);
                this.response=response.data
            } catch(e){
                this.errormsg = e.response.data
            }
            if (!this.errormsg){
                try{
                    let taxonomy = await this.$axios.get(`/taxonomy?search=${this.response.ScientificName}`);
                    this.taxonomy=taxonomy.data
                } catch(e){
                    this.errormsg = e.response.data
                }
            }

            this.dataArray = this.typeof == "proteins" ? this.response.Proteins : this.response.Nucleotides
            let datalink = this.typeof == 'nucleotides' ? '/proteins': '/nucleotides'
            this.otherRouterLink = '/organism/' + this.taxId + datalink
        },
        async findLocus(locus){
            this.loading = true
            if (this.nucleo){
                this.close(this.nucleo.GBSeqLocus)
            }
            try{
                let nucleo = await this.$axios.get(`organism/${this.taxId}/${this.typeof}/${locus}`);
                this.nucleo = this.typeof == "proteins" ? nucleo.data.Proteins[0] : nucleo.data.Nucleotides[0]
            } catch(e){
                this.errormsg = e.response.data
            }
            document.getElementById("popup-"+locus).style.display = "flex"
            this.loading = false
        },
        close(locus){
            document.getElementById("popup-"+locus).style.display = "none"
            this.nucleo = null
        },
        clip(name){
            navigator.clipboard.writeText(name);
            var popup = document.getElementById("myPopup-clip");
            popup.style.visibility = "visible"
            setTimeout(() => popup.style.visibility = "hidden" , 1200);
        }
    },
    mounted(){
        this.taxId = this.$route.params.taxid;
        this.typeof = this.$route.path.split("/").pop();
        this.refresh()
    },
    created(){
        this.$watch(
            () => this.$route.params,
            (toParams) => {
                this.taxId = this.$route.params.taxid;
                this.typeof = this.$route.path.split("/").pop();
                this.refresh()
            })
    },
}


</script>

<template>
<!--TABELLA -->
    <div class="TitlePage">
        {{this.$route.path.split('/').pop() == 'nucleotides' ? 'Nucleotide': 'Protein'}} 
    </div>
    <div class="infoPage">
        of {{ this.response.ScientificName }}
    </div>
    <div class="otherLink">
        Looking for {{this.$route.path.split('/').pop() == 'nucleotides' ? 'proteins': 'nucleotides'}} ? Go to <RouterLink v-if="otherRouterLink" :to="otherRouterLink">{{ this.response.ScientificName }} {{this.$route.path.split('/').pop() == 'nucleotides' ? '/proteins': '/nucleotides'}} </RouterLink>
    </div>
    <ErrorMsg v-if="errormsg" :msg="errormsg"></ErrorMsg>
    <div v-if="this.dataArray">
        <br/>
        <LoadingSpinner :loading="this.loading"/>
        <div style="display: flex;">

            <div class="title-left">
                Locus Id
            </div>
            <div class="title-right">
                Taxonomy
            </div>
        </div>
        <div class="block-organism">
            <div class="left">
                <ul>
                    <div v-for="item in this.dataArray">
                        <div @click="findLocus(item.GBSeq_locus)" class="clickable">
                            {{ item.GBSeq_locus }}
                        </div>
                        <div>
                            <div class="popup" :id='"popup-"+item.GBSeq_locus'>
                                <div v-if="this.nucleo" class="popup-content">
                                    <NucleoComp :data="this.nucleo"/>
                                    <button type="button" class="popupclose" @click="close(item.GBSeq_locus)">Close</button>
                                </div>
                            </div>
                            
                        </div>
                    </div>
                </ul>
            </div>
            <div class="right" v-if="taxonomy">
                <h2>
                    Organism Name: {{ taxonomy.ScientificName }}
                </h2>
                <h2>
                    Taxonomy Id: {{ taxonomy.TaxID }}
                </h2>
                <h2>
                    Rank: {{ taxonomy.Rank }}
                </h2>
                <h2>
                    Division: {{ taxonomy.Division }}
                </h2>
                <h2>
                    Lineage:
                    <ul v-for="lin in taxonomy.LineageEx">
                        <h3>
                            -----------------------------
                        </h3>

                        <RouterLink :to="'/taxonomy/'+lin.TaxID">
                            <h3>
                                Taxonomy Id: {{ lin.TaxID }}
                            </h3>
                            <h3>
                                Scientific Name: {{ lin.ScientificName }}
                            </h3>
                            <h3>
                                Rank: {{ lin.Rank }}
                            </h3>
                        </RouterLink>
                        <h3>
                            -----------------------------
                        </h3>
                    </ul>
                </h2>
                
            </div>
        </div>
        <div v-if="this.response.Products && this.response.Products.length > 0">
            <div class="title-down">
                Products and occurences
            </div>
            <div style="display: flex;">
                <div class="product">
                    <div v-for="prod in this.response.Products" style="display: flex;">
                        <div class="copyToken clickable" @click="clip(prod.ProductName)">
                            <h2>
                                {{ prod.ProductName }}
                            </h2>
                        </div>
                        <h3 style="margin-left: 5%;">
                            {{ prod.QtyProduct }}
                        </h3>
                    </div>
                </div>
                <div class="popup-clip">
                    <span class="popuptext-clip" id="myPopup-clip">Product Copied to ClipBoard!</span>
                </div>
            </div>
        </div>
        <div v-if="this.response.Country && this.response.Country.length > 0">
            <div class="title-down">
                Country
            </div>
            <div style="display: flex;">
                <div class="product">
                    <div v-for="prod in this.response.Country" style="display: flex;">
                        <div class="clickable">
                            <h2>
                                <a :href="'https://www.google.com/maps/place/'+prod.CountryName+'/'" target=”_blank”>
                                    {{ prod.CountryName }}
                                </a>
                            </h2>
                        </div>
                    </div>
                </div>
                <div class="popup-clip">
                    <span class="popuptext-clip" id="myPopup-clip">Product Copied to ClipBoard!</span>
                </div>
            </div>
        </div>
    </div>
    

</template>

<style>
.TitlePage {
    font-size: 38px;
    color: orange;
    max-width: 100px;
}
.TitlePage:hover, .TitlePage:focus, .TitlePage:active { 
    /* color:#ffcc00;  */
    text-shadow: 0 2px 3px #000;
    
}

.infoPage{
    padding-top: 5px;
    font-size: 25px;
}

.otherLink{
    padding-top: 10px;
    font-size: 15px;
}

.clickable:hover, .clickable:focus, .clickable:active { color:#ffcc00; }

.block-organism {
    display: flex;
}
.title-down{
    bottom: -45px;
    margin-top: 5%;
    font-weight: bold;
}
.title-left{
    margin-right: auto;
    font-weight: bold;
}
.title-right{
    margin-left: auto;
    font-weight: bold;
    margin-right: 33%;
}

.block-organism .left{
    background: rgba(0, 0, 0, 0.6);
    /* overflow: auto;
    max-height: 100%; */
    width: 20%;
    max-height: 500px;
    overflow-y: auto;
    justify-content: center;
    align-items: center;
    border: 2px solid black;
}
.block-organism .right {
    width: 40%;
    margin-left: auto;
    max-height: 500px;
    overflow-y: auto;
    background: rgba(0, 0, 0, 0.6);
    /* overflow: auto;
    max-height: 100%; */
    justify-content: center;
    align-items: center;
    border: 2px solid black;
}
.product {
    width: 50%;
    max-height: 500px;
    min-height: 50px;
    justify-content: left;
    align-items: center;
    overflow-y: auto;
    background: rgba(0, 0, 0, 0.6);
    border: 2px solid black;

}

.popup{
    background: rgba(0, 0, 0, 0.6);
    position:absolute;
    /* overflow: auto;
    max-height: 100%; */
    top: 10%;
    left: 18%;
    right: 5%;
    /* right: 5%; */
    display: none;
    justify-content: center;
    align-items: center;
    border: 2px solid black;
}

.popup-content{
    overflow-y: auto;
    word-wrap: break-word;
    overflow-x:visible;
    top: 20%;
    bottom: 20%;
    right: 10%;
    left: 10%;
    width: 80%;
    color: #ffffff;
    size: 10px;
    box-shadow: 0 0 5px #0099FF;
    background:rgba(0, 0, 0, 0.815);
    padding: 2px;
    border-radius: 5px;
    position:fixed;
}


.popupclose{
  position: relative;
  bottom: 0px;
  right: 0px;
  border-radius: 5px;
  background-color: #003399;
}


.popup-close{
    position: absolute;
    top: -15px;
    right: -15px;
    background: #fff;
    height: 20px;
    width: 20px;
    border-radius: 50%;
    box-shadow: 6px 6px 29px -4px rgba(0,0,0,0.75);
    cursor: pointer;
}
.popup-clip {
  position: relative;
  display:block;
  cursor: pointer;
  -webkit-user-select: none;
  -moz-user-select: none;
  -ms-user-select: none;
  user-select: none;
}

/* The actual popup */
.popup-clip .popuptext-clip {
  visibility: hidden;
  width: 160px;
  background-color: #555;
  color: #fff;
  text-align: center;
  top: -67px;
  left: -90px;
  border-radius: 6px;
  padding: 8px 0;
  position: absolute;
  z-index: 1;
}

/* Popup arrow */
.popup-clip .popuptext-clip::after {
  content: "";
  position: absolute;
  top: 100%;
  left: 15%;
  margin-left: -5px;
  border-width: 5px;
  border-style: solid;
  border-color: #555 transparent transparent transparent;
}

/* Toggle this class - hide and show the popup */
.popup-clip .show-clip {
  visibility: visible;
  -webkit-animation: fadeIn 1s;
  animation: fadeIn 1s;
}

/* Add animation (fade in the popup) */
@-webkit-keyframes fadeIn {
  from {opacity: 0;} 
  to {opacity: 1;}
}

@keyframes fadeIn {
  from {opacity: 0;}
  to {opacity:1 ;}
}
</style>
