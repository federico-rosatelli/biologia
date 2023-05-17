package structures

type Session struct {
	Id          string
	IdSession   string
	TimeSession int64
}

type Credentials struct {
	Email    string `json:"email"`
	Password string `json:"password"`
}

type User struct {
	Id            string
	Username      string
	LastLoginTime int64
	UserInfo      struct {
		IpClient        string
		UserAgentString string
		UserAgent       UserAgent
		TimeZoneClient  int64
	}
}

type UserAgent struct {
	UserAgentString string `json:"userAgentString"`
	Name            string `json:"name"`
	Type            string `json:"type"`
	Version         string `json:"version"`
	VersionMajor    string `json:"versionMajor"`
	Device          struct {
		Name  string      `json:"name"`
		Type  string      `json:"type"`
		Brand string      `json:"brand"`
		CPU   interface{} `json:"CPU"`
	} `json:"device"`
	Engine struct {
		Name         string `json:"name"`
		Type         string `json:"type"`
		Version      string `json:"version"`
		VersionMajor string `json:"versionMajor"`
	} `json:"engine"`
	OperatingSystem struct {
		Name         string `json:"name"`
		Type         string `json:"type"`
		Version      string `json:"version"`
		VersionMajor string `json:"versionMajor"`
	} `json:"operatingSystem"`
}
