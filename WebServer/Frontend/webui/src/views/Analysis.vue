<script>
import LoadingSpinner from '../components/LoadingSpinner.vue';


export default {
    data: function () {
        return {
            errormsg: null,
            loading: false,
            bioProj: null,
            bioProjId: this.$axios.params.bioprojid,

        };
    },
    methods: {
        async mounted(){
            try{
                let response = await this.$axios.get(`/sra/` + this.bioProjId);
                this.bioProj = response.data
            } catch(e){
                this.errormsg = e.response.data
            }
        },

    },
    components: { LoadingSpinner }
}

</script>


<template>
    <LoadingSpinner :loading="this.loading"/>
    <ErrorMsg v-if="errormsg" :msg="errormsg"></ErrorMsg>
    <div v-if="bioProj.length > 0">
        <div>
            <h1>{{ bioProj.TaxId }} Analysis</h1>
        </div>

        <div class="legend top-right">
            <span class="pink">pink: case</span>
            <span class="blue">blue: control</span>
        </div>

        <div>
            <ul>
                <div class="image">Heat Map:
                    <img :src="'/sra/' + bioProj.TaxId + '/heatmap'">
                </div>

                <div class="image">Volcano Plot:
                    <img :src="'/sra/' + bioProj.TaxId + '/volcanoplot'">               
                </div>
            </ul>
        </div>
    </div>

</template>


<style>
.legend {
position: fixed;
top: 20px;
right: 20px;
background-color: #ffffff;
padding: 10px;
border: 1px solid #ccc;
border-radius: 4px;
}

.legend span {
display: block;
margin-bottom: 5px;
}

.pink {
color: pink;
}

.blue {
color: blue;
}
</style>

