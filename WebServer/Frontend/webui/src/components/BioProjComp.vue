<script>
export default {
	props: ['data'],
    methods:{
    }
}
</script>

<template>
    <!-- <div>BioProjects of {{ scientificName }}</div> -->
    <tr>
        BioProject
    </tr>
    <ul>
        <tr v-for="proj in data">
            <div class="origin">
                BioProject:
            </div>
            <div class="clickable">
                {{ proj.BioProjectId }}
            </div>
            <td v-for="sample in proj.BioSamples">
                <div class="origin">
                    BioSample: 
                </div>
                <div class="clickable">
                    {{ sample.BioSampleId }}
                </div>
                <ul>
                    <tr>
                        <td v-for="exp in sample.Experiments">
                            <div class="origin">
                                Experiment: 
                            </div>
                            <div class="clickable">
                                {{ exp.ExperimentId }}
                            </div>
                            <ul v-for="run in exp.Runs">
                                <tr>
                                    <td>
                                        <div class="origin">
                                            Run:
                                        </div>
                                        <div class="clickable">
                                            {{ run }}
                                        </div>
                                    </td>
                                </tr>
                            </ul>
                        </td>
                    </tr>
                </ul>
            </td>
        </tr>
    </ul>
    <div v-if="data[0].ImagesPath && data[0].ImagesPath.length > 0" class="immagini">
        <div v-for="img in data[0].ImagesPath">
            <img :src="img" alt="Prova">
        </div>
    </div>
</template>

<style>
.origin{
    color: rgb(172, 172, 172);
}
</style>