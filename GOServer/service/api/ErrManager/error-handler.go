package ErrManager

func ErrorHandler(errorType int) int {
	var errInt int
	switch errorType {
	case 0:
		errInt = 500
	case 1:
		errInt = 401
	case 2:
		errInt = 400
	case 3:
		errInt = 403
	case 4:
		errInt = 404
	case 5:
		errInt = 511
	default:
		errInt = 400
	}
	return errInt
}
