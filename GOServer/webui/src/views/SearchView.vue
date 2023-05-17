<script>
export default {
	data: function() {
		return {
			errormsg: null,
			loading: false,
			some_data: null,
			username: null,
		}
	},
	methods: {
        async searchUser(){
            this.loading = true;
			this.errormsg = null;
			try {
				let response = await this.$axios.get(`/search?query=${this.username}`);
				let response1 = await this.$axios.get(`/organism/${this.taxid}/nucleotides`);
				this.some_data = response.data;
			} catch (e) {
				this.errormsg = e.toString();
			}
        },
	},
}
</script>

<template>
	<div>
		<div
			class="d-flex justify-content-between flex-wrap flex-md-nowrap align-items-center pt-3 pb-2 mb-3 border-bottom">
			<h1 class="h2">Search Profiles</h1>
		</div>
        <div class="login-form">
            <input type="text" placeholder="Username" v-model="username" >
			<input type="submit" value="Search" @click="searchUser">
        </div>
    
        <table id="table">
            <tr>
                <td>Username</td>
            </tr>
            <tr v-for="item in this.some_data" :key="item">
                <h5>{{ this.some_data }}</h5>
            </tr>
        </table>
        

		<ErrorMsg v-if="errormsg" :msg="errormsg"></ErrorMsg>
	</div>
</template>

<style>
.row{
 display: flex;
}
#table {
 font-family: Arial, Helvetica, sans-serif;
 border-collapse: collapse;
 width: 100%;
}
#table td, #table th {
 border: 1px solid #ddd;
 padding: 20px;
}
#table .title{
  font-size: 25px;
  background-color: #999999;
  text-align: center;
}
#table tr:nth-child(even){background-color: #f2f2f2;}
#table tr:hover {background-color: #ddd;}
#table th {
 padding-top: 12px;
 padding-bottom: 12px;
 text-align: left;
 background-color: #04AA6D;
 color: white;
}

#table a{
  color: #000;
}
</style>
