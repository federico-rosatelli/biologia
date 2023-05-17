package ErrManager

type StatusError struct {
	message   string
	typeError int
}

type Errors interface {
	Error() string
	Type() int
}

func NewError(errMessage string, errType int) Errors {
	return &StatusError{
		message:   errMessage,
		typeError: errType,
	}
}

func (e *StatusError) Error() string {
	return e.message
}

func (e *StatusError) Type() int {
	return e.typeError
}
