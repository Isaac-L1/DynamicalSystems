function send_error(dest, task, e, console)

    mail = 'ielhome253@gmail.com';
    password = 'Bananas253';
    setpref('Internet','SMTP_Server','smtp.gmail.com');
    setpref('Internet','E_mail',mail);
    setpref('Internet','SMTP_Username',mail);
    setpref('Internet','SMTP_Password',password);
    props = java.lang.System.getProperties;
    props.setProperty('mail.smtp.auth','true');
    props.setProperty('mail.smtp.starttls.enable','true');
    %props.setProperty('mail.smtp.socketFactory.class', 'javax.net.ssl.SSLSocketFactory');
    props.setProperty('mail.smtp.socketFactory.port','465');
    sendmail(dest, 'An Error Occured', ['Your task: ' task ' has ran into an error. \n The errorID is: ' e.identifier '. \n The error message is: ' e.message '. \n The console log is attached.'], console);

end