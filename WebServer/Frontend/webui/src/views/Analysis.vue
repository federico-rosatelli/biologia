<script>
import LoadingSpinner from '../components/LoadingSpinner.vue';


export default {
    data: function () {
        return {
            errormsg: null,
            loading: false,
            organisms: null,

        };
    },
    methods: {
        async mounted(){
            try{
                let response = await this.$axios.get(`/analysis`);
                this.organismsWG = response.data
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
    <div v-if="response.length > 0">
        <div>
            <h1>{{ response.length }} Results</h1>
        </div>

        <div>
            <ul>
                <li v-for="organism in this.organism" :key="organism.TaxId">

                    <strong>{{ organism.ScientificName }}</strong>

                    <h2>CSV File:</h2>
                        <ul>
                            <li v-for="(values, fieldName) in organism.CsvFile" :key="fieldName">
                                <strong>{{ fieldName }}:</strong>
                                <ul>
                                    <li v-for="value in values" :key="value">
                                        {{ value }}
                                    </li>
                                </ul>
                            </li>
                        </ul>
                        <h2>Images:</h2>
                            <ul>
                                <div class="image">
                                    <img :src="'/images/analysis' + organism.TaxId">
                                </div>
                            </ul>
                </li>
            </ul>
        </div>
    </div>

</template>

<style>
strong {
  font-weight: bold;
}
</style>